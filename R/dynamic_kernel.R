# Stationary Dynamic MCMSEM -------------------------------------------------

.normalize_kernel <- function(kernel) {
  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel)) {
    stop(
      paste0(
        "`kernel` must be exactly one of \"contemporaneous\" or \"dynamic\", ",
        "or the supported alias \"static\"."
      ),
      call. = FALSE
    )
  }
  if (identical(kernel, "static")) {
    return("contemporaneous")
  }
  allowed <- c("contemporaneous", "dynamic")
  if (!(kernel %in% allowed)) {
    stop(
      "Invalid `kernel`: ", encodeString(kernel, quote = "\""),
      paste0(
        ". Use exactly \"contemporaneous\" or \"dynamic\", or the alias ",
        "\"static\"; partial matching is not supported."
      ),
      call. = FALSE
    )
  }
  kernel
}

.model_kernel <- function(model) {
  kernel <- tryCatch(model$meta_data$kernel, error = function(e) NULL)
  if (is.null(kernel) || length(kernel) != 1L || is.na(kernel)) {
    return("contemporaneous")
  }
  .normalize_kernel(kernel)
}

.result_kernel <- function(result) {
  kernel <- tryCatch(result$kernel, error = function(e) NULL)
  if (is.null(kernel) || length(kernel) != 1L || is.na(kernel)) {
    kernel <- tryCatch(result$info$kernel, error = function(e) NULL)
  }
  if (is.null(kernel) || length(kernel) != 1L || is.na(kernel)) {
    return(.model_kernel(result$model))
  }
  .normalize_kernel(kernel)
}

.kernel_label <- function(kernel) {
  kernel <- .normalize_kernel(kernel)
  if (identical(kernel, "dynamic")) {
    "Stationary Dynamic MCMSEM"
  } else {
    "Contemporaneous Structural MCMSEM"
  }
}

.dynamic_parameter_labels <- function(variable_names) {
  p <- length(variable_names)
  out <- matrix("", p, p, dimnames = list(variable_names, variable_names))
  for (row in seq_len(p)) {
    for (col in seq_len(p)) {
      out[row, col] <- if (row == col) {
        paste0("phi_", variable_names[row])
      } else {
        paste0(variable_names[col], "_lag_to_", variable_names[row])
      }
    }
  }
  out
}

.dynamic_cholesky_labels <- function(variable_names) {
  p <- length(variable_names)
  out <- matrix("0", p, p, dimnames = list(variable_names, variable_names))
  for (row in seq_len(p)) {
    for (col in seq_len(row)) {
      out[row, col] <- if (row == col) {
        paste0("log_sd_G_", variable_names[row])
      } else {
        paste0("chol_G_", variable_names[row], "_", variable_names[col])
      }
    }
  }
  out
}

.normalize_dynamic_residual_family <- function(residual_family = NULL,
                                               gaussian_residual = TRUE,
                                               gaussian_residual_missing = FALSE) {
  if (is.null(residual_family)) {
    return(if (isTRUE(gaussian_residual)) "gaussian" else "none")
  }
  if (!is.character(residual_family) || length(residual_family) != 1L ||
      is.na(residual_family)) {
    stop(
      "`residual_family` must be exactly one of \"none\", \"gaussian\", or \"common_gamma\".",
      call. = FALSE
    )
  }
  if (identical(residual_family, "gamma")) residual_family <- "common_gamma"
  allowed <- c("none", "gaussian", "common_gamma")
  if (!(residual_family %in% allowed)) {
    stop(
      "Invalid `residual_family`: ", encodeString(residual_family, quote = "\""),
      ". Use exactly \"none\", \"gaussian\", or \"common_gamma\".",
      call. = FALSE
    )
  }
  if (!isTRUE(gaussian_residual_missing)) {
    if (isTRUE(gaussian_residual) && !identical(residual_family, "gaussian")) {
      stop(
        "`gaussian_residual = TRUE` conflicts with `residual_family = ",
        encodeString(residual_family, quote = "\""), ".",
        call. = FALSE
      )
    }
    if (!isTRUE(gaussian_residual) && identical(residual_family, "gaussian")) {
      stop(
        "`gaussian_residual = FALSE` conflicts with `residual_family = \"gaussian\"`.",
        call. = FALSE
      )
    }
  }
  residual_family
}

.dynamic_residual_family <- function(model) {
  family <- tryCatch(model$meta_data$residual_family, error = function(e) NULL)
  if (is.null(family) || length(family) != 1L || is.na(family)) {
    return(if (isTRUE(model$meta_data$gaussian_residual)) "gaussian" else "none")
  }
  .normalize_dynamic_residual_family(
    family, gaussian_residual = identical(family, "gaussian"),
    gaussian_residual_missing = TRUE
  )
}

.dynamic_gamma_loading_labels <- function(variable_names) {
  matrix(
    paste0("loading_Gamma_", variable_names), length(variable_names), 1L,
    dimnames = list(variable_names, "common_gamma")
  )
}

.nearest_spd_start <- function(x, diagonal_floor = 1e-5) {
  # Used only to obtain an admissible optimization start. The fitted covariance
  # is PSD by construction and is never post-hoc clipped.
  x <- (as.matrix(x) + t(as.matrix(x))) / 2
  eig <- eigen(x, symmetric = TRUE)
  values <- pmax(eig$values, diagonal_floor)
  out <- eig$vectors %*% (values * t(eig$vectors))
  diag(out) <- pmax(diag(out), diagonal_floor)
  (out + t(out)) / 2
}

.covariance_to_cholesky_parameters <- function(Psi, diagonal_floor = 1e-8) {
  Psi <- .nearest_spd_start(Psi, diagonal_floor = diagonal_floor)
  L <- t(chol(Psi))
  diag(L) <- log(pmax(diag(L), diagonal_floor))
  L[upper.tri(L)] <- 0
  L
}

.dynamic_model <- function(data, residual_family = "gaussian") {
  p <- data$meta_data$ncol
  variable_names <- data$meta_data$colnames
  residual_family <- .normalize_dynamic_residual_family(
    residual_family,
    gaussian_residual = identical(residual_family, "gaussian"),
    gaussian_residual_missing = TRUE
  )
  gaussian_residual <- identical(residual_family, "gaussian")
  common_gamma <- identical(residual_family, "common_gamma")

  B_names <- .dynamic_parameter_labels(variable_names)
  B_values <- matrix(0, p, p, dimnames = dimnames(B_names))
  diag(B_values) <- 0.2

  tau_names <- matrix(
    paste0("tau_", variable_names), p, 1,
    dimnames = list(variable_names, "third_cumulant")
  )
  tau_values <- matrix(0, p, 1, dimnames = dimnames(tau_names))
  for (i in seq_len(p)) {
    idx <- .nd_to_2d_idx(p, i, i, i)
    tau_values[i, 1] <- data$M3[idx$x, idx$y]
  }

  kappa_names <- matrix(
    paste0("kappa_", variable_names), p, 1,
    dimnames = list(variable_names, "fourth_cumulant")
  )
  kappa_values <- matrix(0, p, 1, dimnames = dimnames(kappa_names))
  for (i in seq_len(p)) {
    idx <- .nd_to_2d_idx(p, i, i, i, i)
    kappa_values[i, 1] <- data$M4[idx$x, idx$y] - 3 * data$M2[i, i]^2
  }
  kappa_values[, 1] <- pmin(pmax(kappa_values[, 1], -1.5), 20)

  L_names <- matrix("0", p, p, dimnames = list(variable_names, variable_names))
  L_values <- matrix(0, p, p, dimnames = dimnames(L_names))
  if (isTRUE(gaussian_residual)) {
    L_names <- .dynamic_cholesky_labels(variable_names)
    Psi_start <- diag(pmax(diag(data$M2) * 0.10, 1e-4), p)
    L_values <- .covariance_to_cholesky_parameters(Psi_start)
  }

  gamma_loading_names <- matrix(
    "0", p, 1L, dimnames = list(variable_names, "common_gamma")
  )
  gamma_loading_values <- matrix(
    0, p, 1L, dimnames = dimnames(gamma_loading_names)
  )
  gamma_shape_names <- matrix(
    "0", 1L, 1L, dimnames = list("common_gamma", "shape")
  )
  gamma_shape_values <- matrix(
    0, 1L, 1L, dimnames = dimnames(gamma_shape_names)
  )
  if (isTRUE(common_gamma)) {
    gamma_loading_names <- .dynamic_gamma_loading_labels(variable_names)
    loading_magnitude <- sqrt(pmax(diag(data$M2) * 0.10, 1e-4))
    loading_sign <- rep(1, p)
    if (p > 1L) {
      loading_sign[-1L] <- sign(data$M2[-1L, 1L])
      loading_sign[loading_sign == 0] <- 1
    }
    gamma_loading_values[, 1L] <- loading_magnitude * loading_sign
    gamma_shape_names[1L, 1L] <- "shape_Gamma"
    gamma_shape_values[1L, 1L] <- 4
  }

  named_matrices <- list(
    B = B_names,
    Tau = tau_names,
    Kappa = kappa_names,
    L_G = L_names
  )
  num_matrices <- list(
    B = B_values,
    Tau = tau_values,
    Kappa = kappa_values,
    L_G = L_values
  )
  if (isTRUE(common_gamma)) {
    named_matrices$Lambda_Gamma <- gamma_loading_names
    named_matrices$Shape_Gamma <- gamma_shape_names
    num_matrices$Lambda_Gamma <- gamma_loading_values
    num_matrices$Shape_Gamma <- gamma_shape_values
  }
  named_matrices$D2 <- diag(as.character(1), p)
  named_matrices$D2[named_matrices$D2 == "0"] <- "0"
  num_matrices$D2 <- diag(1, p)

  par_names <- unlist(lapply(named_matrices, function(x) {
    vals <- as.vector(x)
    vals[is.na(suppressWarnings(as.numeric(vals)))]
  }), use.names = FALSE)
  par_names <- unique(par_names)

  lower <- upper <- numeric(length(par_names))
  names(lower) <- names(upper) <- par_names
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      nm <- B_names[i, j]
      lower[nm] <- if (i == j) 0 else -0.98
      upper[nm] <- 0.98
    }
    lower[tau_names[i, 1]] <- -100
    upper[tau_names[i, 1]] <- 100
    lower[kappa_names[i, 1]] <- -1.99
    upper[kappa_names[i, 1]] <- 100
  }
  if (isTRUE(gaussian_residual)) {
    for (i in seq_len(p)) {
      for (j in seq_len(i)) {
        nm <- L_names[i, j]
        if (i == j) {
          lower[nm] <- log(1e-6)
          upper[nm] <- log(1e2)
        } else {
          lower[nm] <- -100
          upper[nm] <- 100
        }
      }
    }
  }
  if (isTRUE(common_gamma)) {
    for (nm in as.vector(gamma_loading_names)) {
      lower[nm] <- -100
      upper[nm] <- 100
    }
    lower[gamma_shape_names[1L, 1L]] <- 1e-6
    upper[gamma_shape_names[1L, 1L]] <- Inf
  }
  bounds <- as.data.frame(rbind(L = lower, U = upper), check.names = FALSE)

  start_values <- c(
    as.vector(B_values), as.vector(tau_values), as.vector(kappa_values),
    if (isTRUE(gaussian_residual)) {
      as.vector(L_values)[is.na(suppressWarnings(as.numeric(as.vector(L_names))))]
    } else numeric(),
    if (isTRUE(common_gamma)) {
      c(as.vector(gamma_loading_values), as.vector(gamma_shape_values))
    } else numeric()
  )
  # The model parser defines the authoritative ordering, so these are replaced
  # by parse() during construction. They are supplied only as safe defaults.
  start_values <- data.frame(matrix(start_values, nrow = 1), check.names = FALSE)

  bound_defaults <- list(
    L = list(phi = 0, B = -0.98, tau = -100, kappa = -1.99,
             log = log(1e-6), chol = -100, loading = -100,
             shape = 1e-6),
    U = list(phi = 0.98, B = 0.98, tau = 100, kappa = 100,
             log = log(1e2), chol = 100, loading = 100,
             shape = Inf)
  )

  model <- mcmmodelclass(
    named_matrices = named_matrices,
    num_matrices = num_matrices,
    start_values = mcmstartvaluesclass(start_values),
    bounds = bounds,
    meta_data = list(
      n_obs = data$meta_data$N,
      n_phenotypes = p,
      n_latent = 0,
      bound_defaults = bound_defaults,
      bound_default = bound_defaults,
      weighted = data$meta_data$weighted,
      data_was_scaled = data$meta_data$data_was_scaled,
      scale_data = data$meta_data$scale_data,
      original_colnames = variable_names,
      latent_names = character(),
      kernel = "dynamic",
      gaussian_residual = isTRUE(gaussian_residual),
      residual_family = residual_family,
      innovation_variances = rep(1, p)
    )
  )
  if (isTRUE(common_gamma)) {
    model <- MCMparameter(
      model, gamma_shape_names[1L, 1L], "free", start = 4,
      transform = "positive", lower = 1e-6, overwrite = TRUE
    )
  }
  model
}

.dynamic_kron_power <- function(A, order) {
  if (!is.numeric(order) || length(order) != 1L || order < 1 || order != as.integer(order)) {
    stop("`order` must be a positive integer.", call. = FALSE)
  }
  out <- A
  if (order > 1L) {
    for (i in 2:order) out <- kronecker(out, A)
  }
  out
}

.dynamic_spectral_radius <- function(B) {
  max(Mod(eigen(B, only.values = TRUE)$values))
}

.dynamic_cumulant_propagation <- function(B, D, order) {
  B <- as.matrix(B)
  p <- nrow(B)
  if (ncol(B) != p) stop("`B` must be square.", call. = FALSE)
  if (.dynamic_spectral_radius(B) >= 1) {
    stop("`B` must be stationary (spectral radius < 1).", call. = FALSE)
  }
  D_vec <- as.vector(D)
  expected <- p^order
  if (length(D_vec) != expected) {
    stop("`D` must contain p^order entries.", call. = FALSE)
  }
  solve(diag(expected) - .dynamic_kron_power(B, order), D_vec)
}

.dynamic_array_to_mcm <- function(A) {
  dims <- dim(A)
  p <- dims[1]
  order <- length(dims)
  if (is.null(dims) || !all(dims == p) || order < 2L) {
    stop("`A` must be an order-2-or-higher cubical array.", call. = FALSE)
  }
  out <- matrix(0, p, p^(order - 1L))
  grid <- as.matrix(expand.grid(rep(list(seq_len(p)), order)))
  for (i in seq_len(nrow(grid))) {
    coord <- as.integer(grid[i, ])
    idx <- do.call(.nd_to_2d_idx, c(list(p, coord[1], coord[2]), as.list(coord[-c(1, 2)])))
    out[idx$x, idx$y] <- A[matrix(coord, nrow = 1)]
  }
  out
}

.dynamic_mcm_to_array <- function(M, order) {
  M <- as.matrix(M)
  p <- nrow(M)
  if (ncol(M) != p^(order - 1L)) {
    stop("Moment matrix dimensions do not match `order`.", call. = FALSE)
  }
  out <- array(0, dim = rep(p, order))
  grid <- as.matrix(expand.grid(rep(list(seq_len(p)), order)))
  for (i in seq_len(nrow(grid))) {
    coord <- as.integer(grid[i, ])
    idx <- do.call(.nd_to_2d_idx, c(list(p, coord[1], coord[2]), as.list(coord[-c(1, 2)])))
    out[matrix(coord, nrow = 1)] <- M[idx$x, idx$y]
  }
  out
}

.cumulant4_to_raw4 <- function(K4, Sigma) {
  if (is.matrix(K4)) K4 <- .dynamic_mcm_to_array(K4, 4L)
  Sigma <- as.matrix(Sigma)
  p <- nrow(Sigma)
  out <- K4
  for (i in seq_len(p)) for (j in seq_len(p)) {
    for (k in seq_len(p)) for (l in seq_len(p)) {
      out[i, j, k, l] <- K4[i, j, k, l] +
        Sigma[i, j] * Sigma[k, l] +
        Sigma[i, k] * Sigma[j, l] +
        Sigma[i, l] * Sigma[j, k]
    }
  }
  out
}

.raw4_to_cumulant4 <- function(M4, Sigma) {
  if (is.matrix(M4)) M4 <- .dynamic_mcm_to_array(M4, 4L)
  Sigma <- as.matrix(Sigma)
  p <- nrow(Sigma)
  out <- M4
  for (i in seq_len(p)) for (j in seq_len(p)) {
    for (k in seq_len(p)) for (l in seq_len(p)) {
      out[i, j, k, l] <- M4[i, j, k, l] -
        Sigma[i, j] * Sigma[k, l] -
        Sigma[i, k] * Sigma[j, l] -
        Sigma[i, l] * Sigma[j, k]
    }
  }
  out
}

.dynamic_gaussian_covariance <- function(L_parameters) {
  L_parameters <- as.matrix(L_parameters)
  if (nrow(L_parameters) != ncol(L_parameters)) {
    stop("Cholesky parameter matrix must be square.", call. = FALSE)
  }
  L <- L_parameters
  L[upper.tri(L)] <- 0
  diag(L) <- exp(diag(L))
  L %*% t(L)
}

.dynamic_vector_outer_power <- function(x, order) {
  x <- as.numeric(x)
  if (!length(x) || length(order) != 1L || order < 1L ||
      order != as.integer(order)) {
    stop("A non-empty vector and positive integer `order` are required.",
         call. = FALSE)
  }
  out <- x
  if (order > 1L) {
    for (i in 2:order) out <- kronecker(out, x)
  }
  array(out, dim = rep(length(x), order))
}

.dynamic_residual_cumulants <- function(model) {
  p <- model$meta_data$n_phenotypes
  family <- .dynamic_residual_family(model)
  zero2 <- matrix(0, p, p)
  zero3 <- array(0, rep(p, 3L))
  zero4 <- array(0, rep(p, 4L))
  if (identical(family, "none")) {
    return(list(
      family = family, M2 = zero2, C3 = zero3, K4 = zero4,
      loadings = numeric(), shape = NA_real_, skewness = 0,
      excess_kurtosis = 0, L_G = zero2, Psi_G = zero2
    ))
  }
  if (identical(family, "gaussian")) {
    Psi_G <- .dynamic_gaussian_covariance(model$num_matrices$L_G)
    return(list(
      family = family, M2 = Psi_G, C3 = zero3, K4 = zero4,
      loadings = numeric(), shape = NA_real_, skewness = 0,
      excess_kurtosis = 0, L_G = model$num_matrices$L_G, Psi_G = Psi_G
    ))
  }
  if (!identical(family, "common_gamma")) {
    stop("Unsupported dynamic residual family `", family, "`.", call. = FALSE)
  }
  loadings <- as.numeric(model$num_matrices$Lambda_Gamma)
  shape <- as.numeric(model$num_matrices$Shape_Gamma)[1L]
  if (length(loadings) != p || !is.finite(shape) || shape <= 0) {
    stop("The common-gamma residual requires finite loadings and positive shape.",
         call. = FALSE)
  }
  skewness <- 2 / sqrt(shape)
  excess_kurtosis <- 6 / shape
  list(
    family = family,
    M2 = tcrossprod(loadings),
    C3 = skewness * .dynamic_vector_outer_power(loadings, 3L),
    K4 = excess_kurtosis * .dynamic_vector_outer_power(loadings, 4L),
    loadings = loadings,
    shape = shape,
    skewness = skewness,
    excess_kurtosis = excess_kurtosis,
    L_G = zero2,
    Psi_G = zero2
  )
}

.dynamic_implied_moments_base <- function(model, parameters = NULL,
                                          stationarity_limit = 1) {
  if (!identical(.model_kernel(model), "dynamic")) {
    stop("A dynamic MCMSEM model is required.", call. = FALSE)
  }
  model <- model$copy()
  if (!is.null(parameters)) {
    if (length(parameters) != length(model$param_values)) {
      stop("Parameter vector has the wrong length.", call. = FALSE)
    }
    model$param_values <- as.numeric(parameters)
    model$inverse_parse()
  }
  B <- model$num_matrices$B
  rho <- .dynamic_spectral_radius(B)
  if (!is.finite(rho) || rho >= stationarity_limit) {
    stop("`B` is not admissibly stationary.", call. = FALSE)
  }
  p <- nrow(B)
  tau <- as.numeric(model$num_matrices$Tau)
  kappa <- as.numeric(model$num_matrices$Kappa)
  D2 <- diag(1, p)
  D3 <- array(0, rep(p, 3L))
  D4 <- array(0, rep(p, 4L))
  for (i in seq_len(p)) {
    D3[matrix(rep(i, 3L), nrow = 1)] <- tau[i]
    D4[matrix(rep(i, 4L), nrow = 1)] <- kappa[i]
  }
  C2 <- array(.dynamic_cumulant_propagation(B, D2, 2L), c(p, p))
  dynamic_C3 <- array(.dynamic_cumulant_propagation(B, D3, 3L), rep(p, 3L))
  dynamic_K4 <- array(.dynamic_cumulant_propagation(B, D4, 4L), rep(p, 4L))
  residual <- .dynamic_residual_cumulants(model)
  C3 <- dynamic_C3 + residual$C3
  K4 <- dynamic_K4 + residual$K4
  Sigma <- C2 + residual$M2
  M4 <- .cumulant4_to_raw4(K4, Sigma)
  dimnames(B) <- list(model$meta_data$original_colnames,
                      model$meta_data$original_colnames)
  dimnames(residual$Psi_G) <- list(model$meta_data$original_colnames,
                                   model$meta_data$original_colnames)
  dimnames(residual$M2) <- list(model$meta_data$original_colnames,
                                model$meta_data$original_colnames)
  if (length(residual$loadings)) {
    names(residual$loadings) <- model$meta_data$original_colnames
  }
  list(
    B = B,
    M2 = Sigma,
    M3 = .dynamic_array_to_mcm(C3),
    K4 = .dynamic_array_to_mcm(K4),
    M4 = .dynamic_array_to_mcm(M4),
    within_M2 = C2,
    dynamic_M3 = .dynamic_array_to_mcm(dynamic_C3),
    dynamic_K4 = .dynamic_array_to_mcm(dynamic_K4),
    residual_M2 = residual$M2,
    residual_M3 = .dynamic_array_to_mcm(residual$C3),
    residual_K4 = .dynamic_array_to_mcm(residual$K4),
    residual_family = residual$family,
    residual_loadings = residual$loadings,
    residual_shape = residual$shape,
    residual_skewness = residual$skewness,
    residual_excess_kurtosis = residual$excess_kurtosis,
    Psi_G = residual$Psi_G,
    L_G = residual$L_G,
    D2 = D2,
    D3 = .dynamic_array_to_mcm(D3),
    D4 = .dynamic_array_to_mcm(D4),
    spectral_radius = rho,
    stationary = rho < 1
  )
}

.dynamic_unique_indices <- function(p, order) {
  grid <- as.matrix(expand.grid(rep(list(seq_len(p)), order)))
  grid[apply(grid, 1, function(x) all(diff(x) >= 0)), , drop = FALSE]
}

.dynamic_unique_moments <- function(M, order) {
  M <- as.matrix(M)
  p <- nrow(M)
  grid <- .dynamic_unique_indices(p, order)
  vapply(seq_len(nrow(grid)), function(i) {
    coord <- as.integer(grid[i, ])
    idx <- do.call(.nd_to_2d_idx, c(list(p, coord[1], coord[2]), as.list(coord[-c(1, 2)])))
    M[idx$x, idx$y]
  }, numeric(1))
}

MCMmomentcount <- function(x, use_skewness = TRUE, use_kurtosis = TRUE) {
  p <- if (length(x) == 1L && is.numeric(x)) {
    as.integer(x)
  } else if (inherits(x, "mcmresultclass")) {
    x$model$meta_data$n_phenotypes
  } else if (inherits(x, "mcmmodelclass")) {
    x$meta_data$n_phenotypes
  } else if (inherits(x, "mcmdataclass")) {
    x$meta_data$ncol
  } else {
    ncol(x)
  }
  if (!is.finite(p) || p < 1L) stop("Could not determine a positive number of variables.")
  counts <- c(
    covariance = choose(p + 1L, 2L),
    third = if (use_skewness) choose(p + 2L, 3L) else 0,
    fourth = if (use_kurtosis) choose(p + 3L, 4L) else 0
  )
  c(counts, total = sum(counts))
}

MCMdegreesoffreedom <- function(object, use_skewness = TRUE,
                                use_kurtosis = TRUE) {
  model <- if (inherits(object, "mcmresultclass")) object$model else object
  if (!inherits(model, "mcmmodelclass")) {
    stop("`object` must be an MCM model or result.", call. = FALSE)
  }
  if (inherits(object, "mcmresultclass")) {
    use_skewness <- isTRUE(object$info$use_skewness)
    use_kurtosis <- isTRUE(object$info$use_kurtosis)
  }
  counts <- MCMmomentcount(model, use_skewness, use_kurtosis)
  n_parameters <- length(model$param_values)
  .ensure_parameter_graph(model)
  list(
    moments_by_order = counts[c("covariance", "third", "fourth")],
    n_moments = unname(counts[["total"]]),
    n_parameters = n_parameters,
    n_reported_parameters = nrow(model$parameter_table),
    n_fixed_parameters = sum(model$parameter_table$type == "fixed"),
    n_derived_parameters = sum(model$parameter_table$type == "derived"),
    df = unname(counts[["total"]]) - n_parameters
  )
}

MCMdiagnostics <- function(object, jacobian = FALSE,
                           tolerance = sqrt(.Machine$double.eps)) {
  model <- if (inherits(object, "mcmresultclass")) object$model else object
  if (!inherits(model, "mcmmodelclass")) {
    stop("`object` must be an MCM model or result.", call. = FALSE)
  }
  out <- MCMdegreesoffreedom(object)
  out$kernel <- .model_kernel(model)
  out$observed_moment_count <- out$n_moments
  out$independent_free_parameter_count <- out$n_parameters
  out$nominal_df <- out$df
  derived <- model$parameter_table$type == "derived"
  out$derived_constraints <- model$parameter_table[
    derived, c("name", "expression"), drop = FALSE
  ]
  if (identical(out$kernel, "dynamic")) {
    implied <- .dynamic_implied_moments_base(model)
    out$spectral_radius <- implied$spectral_radius
    out$stationary <- implied$stationary
  }
  if (isTRUE(jacobian)) {
    optimizer_coordinates <- .parameter_optimizer_coordinates(model)
    map <- function(par) {
      x <- MCMimpliedmoments(
        model, par, parameter_scale = "optimizer",
        use_skewness = if (inherits(object, "mcmresultclass")) {
          isTRUE(object$info$use_skewness)
        } else TRUE,
        use_kurtosis = if (inherits(object, "mcmresultclass")) {
          isTRUE(object$info$use_kurtosis)
        } else TRUE
      )
      c(
        .dynamic_unique_moments(x$M2, 2L),
        if (!is.null(x$M3)) .dynamic_unique_moments(x$M3, 3L) else numeric(),
        if (!is.null(x$M4)) .dynamic_unique_moments(x$M4, 4L) else numeric()
      )
    }
    J <- if (length(optimizer_coordinates)) {
      numDeriv::jacobian(map, optimizer_coordinates)
    } else {
      matrix(numeric(), length(map(optimizer_coordinates)), 0L)
    }
    colnames(J) <- model$param_names
    decomposition <- if (ncol(J)) svd(J) else NULL
    singular_values <- if (is.null(decomposition)) numeric() else decomposition$d
    out$jacobian <- J
    out$jacobian_rank <- qr(J, tol = tolerance)$rank
    out$locally_full_column_rank <- out$jacobian_rank == length(model$param_values)
    out$jacobian_singular_values <- singular_values
    out$jacobian_condition <- if (!ncol(J)) {
      NA_real_
    } else if (length(singular_values) && min(singular_values) > 0) {
      max(singular_values) / min(singular_values)
    } else Inf
    threshold <- tolerance * if (length(singular_values)) max(singular_values) else 0
    near <- which(singular_values <= threshold)
    out$near_null_directions <- if (length(near)) {
      directions <- decomposition$v[, near, drop = FALSE]
      rownames(directions) <- model$param_names
      colnames(directions) <- paste0("direction_", seq_len(ncol(directions)))
      directions
    } else matrix(numeric(), nrow = length(model$param_names), ncol = 0L,
                  dimnames = list(model$param_names, character()))
  }
  out
}
