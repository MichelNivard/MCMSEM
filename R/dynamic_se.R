# WLS weighting and asymptotic covariance for the dynamic kernel ------------

.normalize_moment_weighting <- function(moment_weighting) {
  if (!is.character(moment_weighting) || length(moment_weighting) != 1L ||
      is.na(moment_weighting) ||
      !(moment_weighting %in% c("identity", "diagonal", "full"))) {
    stop(
      "`moment_weighting` must be exactly one of \"identity\", \"diagonal\", or \"full\".",
      call. = FALSE
    )
  }
  moment_weighting
}

.normalize_se_correction <- function(se_correction) {
  if (!is.character(se_correction) || length(se_correction) != 1L ||
      is.na(se_correction) ||
      !(se_correction %in% c("auto", "robust", "model_based"))) {
    stop(
      "`se_correction` must be exactly one of \"auto\", \"robust\", or \"model_based\".",
      call. = FALSE
    )
  }
  se_correction
}

.dynamic_moment_vector <- function(moment_list) {
  c(
    .dynamic_unique_moments(moment_list$M2, 2L),
    .dynamic_unique_moments(moment_list$M3, 3L),
    .dynamic_unique_moments(moment_list$M4, 4L)
  )
}

.dynamic_moment_grids <- function(p) {
  list(
    M2 = .dynamic_unique_indices(p, 2L),
    M3 = .dynamic_unique_indices(p, 3L),
    M4 = .dynamic_unique_indices(p, 4L)
  )
}

.dynamic_torch_unique_moments <- function(moment, grid) {
  p <- moment$shape[[1]]
  entries <- lapply(seq_len(nrow(grid)), function(i) {
    coord <- as.integer(grid[i, ])
    idx <- do.call(
      .nd_to_2d_idx,
      c(list(p, coord[1], coord[2]), as.list(coord[-c(1, 2)]))
    )
    moment[idx$x, idx$y]
  })
  torch_stack(entries)
}

.dynamic_torch_moment_vector <- function(moment_list, grids) {
  torch_cat(list(
    .dynamic_torch_unique_moments(moment_list$M2, grids$M2),
    .dynamic_torch_unique_moments(moment_list$M3, grids$M3),
    .dynamic_torch_unique_moments(moment_list$M4, grids$M4)
  ))
}

.dynamic_require_moment_vcov <- function(data) {
  if (!isTRUE(data$SE$computed) || is.null(data$SE$S.m)) {
    stop(
      "Dynamic WLS weighting or asymptotic SEs require a data summary made with `prep_asymptotic_se = TRUE`.",
      call. = FALSE
    )
  }
  if (!isTRUE(data$SE$influence_function_corrected) ||
      !identical(data$SE$representation, "raw_central_moments")) {
    stop(
      paste0(
        "This summary does not contain the mean-corrected covariance of raw central moments required by the dynamic kernel. ",
        "Recreate it with MCMSEM >= 0.27.0 and `prep_asymptotic_se = TRUE`."
      ),
      call. = FALSE
    )
  }
  if (isTRUE(data$meta_data$weighted)) {
    stop(
      "Dynamic WLS weighting and asymptotic SEs do not yet support analysis weights.",
      call. = FALSE
    )
  }
  expected <- unname(MCMmomentcount(data)[["total"]])
  omega <- .dynamic_tensor_to_matrix(data$SE$S.m)
  idx <- data$SE$idx$idx
  if (!is.null(idx)) omega <- omega[idx, idx, drop = FALSE]
  if (length(dim(omega)) != 2L || any(dim(omega) != expected) ||
      any(!is.finite(omega))) {
    stop("The prepared moment covariance has incompatible dimensions or non-finite values.", call. = FALSE)
  }
  (omega + t(omega)) / 2
}

.dynamic_regularized_inverse <- function(x, ridge = 1e-8) {
  x <- (as.matrix(x) + t(as.matrix(x))) / 2
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("`weight_ridge` must be one finite, non-negative number.", call. = FALSE)
  }
  eig <- eigen(x, symmetric = TRUE)
  reference <- max(abs(eig$values), .Machine$double.eps)
  floor_value <- max(ridge * reference, .Machine$double.eps * reference)
  regularized_values <- pmax(eig$values, floor_value)
  regularized <- eig$vectors %*% (regularized_values * t(eig$vectors))
  inverse <- eig$vectors %*% ((1 / regularized_values) * t(eig$vectors))
  list(
    inverse = (inverse + t(inverse)) / 2,
    regularized = (regularized + t(regularized)) / 2,
    eigenvalues = eig$values,
    floor = floor_value,
    condition = max(regularized_values) / min(regularized_values)
  )
}

.dynamic_weight_specification <- function(data, moment_weighting = "identity",
                                          weight_ridge = 1e-8,
                                          require_vcov = FALSE) {
  moment_weighting <- .normalize_moment_weighting(moment_weighting)
  if (length(weight_ridge) != 1L || !is.finite(weight_ridge) ||
      weight_ridge < 0) {
    stop("`weight_ridge` must be one finite, non-negative number.", call. = FALSE)
  }
  n_moments <- unname(MCMmomentcount(data)[["total"]])
  omega <- NULL
  if (isTRUE(require_vcov) || !identical(moment_weighting, "identity")) {
    omega <- .dynamic_require_moment_vcov(data)
  }
  if (identical(moment_weighting, "identity")) {
    W <- diag(n_moments)
    regularization <- NULL
  } else if (identical(moment_weighting, "diagonal")) {
    diagonal <- diag(omega)
    reference <- max(abs(diagonal), .Machine$double.eps)
    floor_value <- max(weight_ridge * reference, .Machine$double.eps * reference)
    W <- diag(1 / pmax(diagonal, floor_value), n_moments)
    regularization <- list(floor = floor_value)
  } else {
    regularization <- .dynamic_regularized_inverse(omega, weight_ridge)
    W <- regularization$inverse
  }
  # A scalar multiple of W has the same minimizer and asymptotic covariance.
  # Normalization keeps optimizer magnitudes comparable across N and choices of W.
  normalization <- mean(diag(W))
  if (!is.finite(normalization) || normalization <= 0) {
    stop("The moment weight matrix is not positive definite.", call. = FALSE)
  }
  list(
    type = moment_weighting,
    W = (W / normalization + t(W / normalization)) / 2,
    omega = omega,
    normalization = normalization,
    regularization = regularization
  )
}

.dynamic_gaussian_vector <- function(model, parameters) {
  reported <- .parameter_values_base(model, parameters, optimizer_scale = TRUE)
  Psi <- .dynamic_implied_moments_base(
    model, unname(reported[model$param_names])
  )$Psi_G
  unlist(lapply(seq_len(nrow(Psi)), function(row) Psi[row, seq_len(row)]),
         use.names = FALSE)
}

.dynamic_asymptotic_se <- function(model, data, weight_spec,
                                   se_correction = "auto",
                                   jacobian_method = "simple",
                                   weight_ridge = 1e-8) {
  se_correction <- .normalize_se_correction(se_correction)
  if (is.null(weight_spec$omega)) {
    weight_spec$omega <- .dynamic_require_moment_vcov(data)
  }
  theta <- .parameter_optimizer_coordinates(model)
  implied <- function(optimizer_coordinates) {
    reported <- .parameter_values_base(
      model, optimizer_coordinates, optimizer_scale = TRUE
    )
    free_values <- unname(reported[model$param_names])
    .dynamic_moment_vector(.dynamic_implied_moments_base(model, free_values))
  }
  if (length(theta) == 0L) {
    all_names <- model$parameter_table$name
    Delta <- matrix(numeric(), length(implied(theta)), 0L,
                    dimnames = list(NULL, character()))
    V <- matrix(0, length(all_names), length(all_names),
                dimnames = list(all_names, all_names))
    correction <- if (identical(se_correction, "auto")) {
      if (identical(weight_spec$type, "full")) "model_based" else "robust"
    } else se_correction
    gaussian_count <- if (isTRUE(model$meta_data$gaussian_residual)) {
      model$meta_data$n_phenotypes * (model$meta_data$n_phenotypes + 1L) / 2L
    } else 0L
    return(list(
      se = stats::setNames(rep(0, length(all_names)), all_names),
      vcov = V, vcov_robust = V, vcov_model_based = V,
      vcov_optimizer = matrix(numeric(), 0L, 0L),
      vcov_optimizer_robust = matrix(numeric(), 0L, 0L),
      vcov_optimizer_model_based = matrix(numeric(), 0L, 0L),
      vcov_free = matrix(numeric(), 0L, 0L),
      parameter_jacobian = matrix(0, length(all_names), 0L,
                                  dimnames = list(all_names, character())),
      jacobian = Delta, jacobian_rank = 0L,
      jacobian_singular_values = numeric(), jacobian_condition = NA_real_,
      information = matrix(numeric(), 0L, 0L), bread_condition = NA_real_,
      correction = correction,
      gaussian_jacobian = matrix(numeric(), gaussian_count, 0L),
      gaussian_vcov = matrix(0, gaussian_count, gaussian_count),
      gaussian_se = rep(0, gaussian_count)
    ))
  }
  Delta <- numDeriv::jacobian(
    func = implied, x = theta, method = jacobian_method
  )
  colnames(Delta) <- model$param_names
  moment_names <- c(
    paste0("M2_", apply(.dynamic_unique_indices(data$meta_data$ncol, 2L), 1, paste, collapse = "_")),
    paste0("M3_", apply(.dynamic_unique_indices(data$meta_data$ncol, 3L), 1, paste, collapse = "_")),
    paste0("M4_", apply(.dynamic_unique_indices(data$meta_data$ncol, 4L), 1, paste, collapse = "_"))
  )
  rownames(Delta) <- moment_names
  singular_values <- svd(Delta, nu = 0, nv = 0)$d
  jacobian_condition <- max(singular_values) / min(singular_values)
  rank <- qr(Delta, tol = sqrt(.Machine$double.eps))$rank
  if (rank < length(theta)) {
    warning(
      "The dynamic moment Jacobian is rank deficient; asymptotic standard errors are not identified at this solution.",
      call. = FALSE
    )
    all_names <- model$parameter_table$name
    V_na <- matrix(NA_real_, length(all_names), length(all_names),
                   dimnames = list(all_names, all_names))
    V_coordinate_na <- matrix(
      NA_real_, length(theta), length(theta),
      dimnames = list(model$param_names, model$param_names)
    )
    return(list(
      se = stats::setNames(rep(NA_real_, length(all_names)), all_names),
      vcov = V_na, vcov_robust = V_na, vcov_model_based = V_na,
      vcov_optimizer = V_coordinate_na, vcov_free = V_coordinate_na,
      vcov_optimizer_robust = V_coordinate_na,
      vcov_optimizer_model_based = V_coordinate_na,
      jacobian = Delta, jacobian_rank = rank,
      jacobian_singular_values = singular_values,
      jacobian_condition = jacobian_condition, bread_condition = Inf,
      correction = NA_character_, gaussian_vcov = matrix(),
      gaussian_se = numeric()
    ))
  }

  W <- weight_spec$W
  omega <- weight_spec$omega
  information <- crossprod(Delta, W %*% Delta)
  bread_spec <- .dynamic_regularized_inverse(information, weight_ridge)
  bread <- bread_spec$inverse
  meat <- crossprod(Delta, W %*% omega %*% W %*% Delta)
  robust <- bread %*% meat %*% bread
  robust <- (robust + t(robust)) / 2
  # For full WLS, the unnormalized W is the regularized inverse of Var(s).
  # Undo the scalar optimizer normalization before applying (Delta' W Delta)^-1.
  model_based <- (bread + t(bread)) / (2 * weight_spec$normalization)
  correction <- if (identical(se_correction, "auto")) {
    if (identical(weight_spec$type, "full")) "model_based" else "robust"
  } else {
    se_correction
  }
  V_optimizer <- if (identical(correction, "robust")) robust else model_based
  dimnames(V_optimizer) <- dimnames(robust) <- dimnames(model_based) <-
    list(model$param_names, model$param_names)
  reported_covariance <- .parameter_covariance_from_optimizer(
    model, V_optimizer, theta, method = jacobian_method
  )
  robust_reported <- .parameter_covariance_from_optimizer(
    model, robust, theta, method = jacobian_method
  )$vcov
  model_based_reported <- .parameter_covariance_from_optimizer(
    model, model_based, theta, method = jacobian_method
  )$vcov
  V <- reported_covariance$vcov
  se <- reported_covariance$se

  if (isTRUE(model$meta_data$gaussian_residual)) {
    gaussian_jacobian <- numDeriv::jacobian(
      func = function(parameters) .dynamic_gaussian_vector(model, parameters),
      x = theta, method = jacobian_method
    )
    gaussian_vcov <- gaussian_jacobian %*% V_optimizer %*% t(gaussian_jacobian)
    gaussian_vcov <- (gaussian_vcov + t(gaussian_vcov)) / 2
    gaussian_se <- sqrt(pmax(diag(gaussian_vcov), 0))
  } else {
    gaussian_jacobian <- matrix(numeric(), 0L, length(theta))
    gaussian_vcov <- matrix(numeric(), 0L, 0L)
    gaussian_se <- numeric()
  }

  list(
    se = se,
    vcov = V,
    vcov_robust = robust_reported,
    vcov_model_based = model_based_reported,
    vcov_optimizer = V_optimizer,
    vcov_optimizer_robust = robust,
    vcov_optimizer_model_based = model_based,
    vcov_free = reported_covariance$free_vcov,
    parameter_jacobian = reported_covariance$jacobian,
    jacobian = Delta,
    jacobian_rank = rank,
    jacobian_singular_values = singular_values,
    jacobian_condition = jacobian_condition,
    information = information,
    bread_condition = bread_spec$condition,
    correction = correction,
    gaussian_jacobian = gaussian_jacobian,
    gaussian_vcov = gaussian_vcov,
    gaussian_se = gaussian_se
  )
}
