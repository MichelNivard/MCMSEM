.dynamic_model_with_parameters <- function(model, parameters) {
  out <- model$copy()
  out$param_values <- as.numeric(parameters)
  out$start_values$set_all(as.numeric(parameters))
  out$inverse_parse()
  out
}

.dynamic_screen_start <- function(model, parameters,
                                  stationarity_limit = 0.995) {
  out <- as.numeric(parameters)
  names(out) <- model$param_names
  lower <- as.numeric(model$bounds["L", model$param_names, drop = TRUE])
  upper <- as.numeric(model$bounds["U", model$param_names, drop = TRUE])
  out <- pmin(pmax(out, lower + 1e-8), upper - 1e-8)
  candidate <- .dynamic_model_with_parameters(model, out)
  B <- candidate$num_matrices$B
  rho <- .dynamic_spectral_radius(B)
  if (!is.finite(rho)) stop("Starting `B` has non-finite eigenvalues.", call. = FALSE)
  if (rho >= stationarity_limit) {
    B <- B * (0.90 * stationarity_limit / rho)
    B_names <- candidate$named_matrices$B
    for (i in seq_along(B_names)) {
      nm <- gsub("-", "", B_names[i])
      if (nm %in% names(out)) out[[nm]] <- B[i]
    }
  }
  unname(out)
}

.dynamic_random_starts <- function(model, n_starts, seed,
                                   stationarity_limit = 0.995) {
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L) {
    stop("`n_starts` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(seed)
    if (exists("torch_manual_seed", mode = "function")) {
      torch_manual_seed(as.integer(seed))
    }
  }
  starts <- vector("list", n_starts)
  starts[[1L]] <- .dynamic_screen_start(
    model, model$param_values, stationarity_limit
  )
  if (n_starts == 1L) return(starts)

  B_names <- as.vector(model$named_matrices$B)
  p <- nrow(model$named_matrices$B)
  tau_names <- as.vector(model$named_matrices$Tau)
  kappa_names <- as.vector(model$named_matrices$Kappa)
  L_names <- as.vector(model$named_matrices$L_G)
  L_names <- L_names[is.na(suppressWarnings(as.numeric(L_names)))]

  for (s in 2:n_starts) {
    start <- model$param_values
    names(start) <- model$param_names
    B <- matrix(0, p, p)
    diag(B) <- stats::runif(p, 0.05, 0.70)
    B[row(B) != col(B)] <- stats::runif(p * (p - 1L), -0.30, 0.30)
    for (i in seq_along(B_names)) start[[B_names[i]]] <- B[i]
    start[tau_names] <- start[tau_names] + stats::rnorm(p, 0, 0.5)
    start[kappa_names] <- start[kappa_names] + stats::rnorm(p, 0, 1)
    if (length(L_names) > 0L) {
      diagonal_L <- grep("^log_sd_G_", L_names, value = TRUE)
      off_diagonal_L <- setdiff(L_names, diagonal_L)
      start[diagonal_L] <- start[diagonal_L] + stats::rnorm(length(diagonal_L), 0, 0.5)
      if (length(off_diagonal_L) > 0L) {
        start[off_diagonal_L] <- start[off_diagonal_L] +
          stats::rnorm(length(off_diagonal_L), 0, 0.25)
      }
    }
    starts[[s]] <- .dynamic_screen_start(
      model, start, stationarity_limit
    )
  }
  starts
}

.dynamic_parameters_admissible <- function(model, parameters, use_bounds,
                                           stationarity_limit,
                                           tolerance = 1e-6) {
  if (any(!is.finite(parameters))) {
    return(list(admissible = FALSE, spectral_radius = NA_real_,
                stationary = FALSE, bounds_ok = FALSE))
  }
  candidate <- .dynamic_model_with_parameters(model, parameters)
  rho <- .dynamic_spectral_radius(candidate$num_matrices$B)
  lower <- as.numeric(candidate$bounds["L", candidate$param_names, drop = TRUE])
  upper <- as.numeric(candidate$bounds["U", candidate$param_names, drop = TRUE])
  bounds_ok <- !isTRUE(use_bounds) || all(parameters >= lower - tolerance &
                                           parameters <= upper + tolerance)
  stationary <- is.finite(rho) && rho < stationarity_limit
  list(
    admissible = stationary && bounds_ok,
    spectral_radius = rho,
    stationary = stationary,
    bounds_ok = bounds_ok
  )
}

.dynamic_tensor_to_matrix <- function(x) {
  as.matrix(torch_tensor(x, device = torch_device("cpu")))
}

.dynamic_tensor_to_numeric <- function(x) {
  as.numeric(torch_tensor(x, device = torch_device("cpu")))
}

.dynamic_transition_table <- function(model, B) {
  p <- nrow(B)
  variables <- model$meta_data$original_colnames
  rows <- vector("list", p^2)
  iter <- 1L
  for (row in seq_len(p)) for (col in seq_len(p)) {
    rows[[iter]] <- data.frame(
      label = model$named_matrices$B[row, col],
      current = variables[row],
      lagged = variables[col],
      type = if (row == col) "autoregressive" else "cross-lagged",
      estimate = B[row, col],
      stringsAsFactors = FALSE
    )
    iter <- iter + 1L
  }
  do.call(rbind, rows)
}

.dynamic_gaussian_table <- function(Psi_G, variable_names,
                                    standard_errors = NULL) {
  p <- nrow(Psi_G)
  rows <- list()
  iter <- 1L
  for (row in seq_len(p)) for (col in seq_len(row)) {
    rows[[iter]] <- data.frame(
      label = if (row == col) {
        paste0("Var(G_", variable_names[row], ")")
      } else {
        paste0("Cov(G_", variable_names[row], ",G_", variable_names[col], ")")
      },
      lhs = variable_names[row],
      rhs = variable_names[col],
      estimate = Psi_G[row, col],
      se = if (is.null(standard_errors)) NA_real_ else standard_errors[iter],
      stringsAsFactors = FALSE
    )
    iter <- iter + 1L
  }
  do.call(rbind, rows)
}

.dynamic_gradient_history_frame <- function(history) {
  if (length(history) == 0L) return(data.frame())
  flattened <- lapply(history, function(snapshot) {
    unlist(snapshot, recursive = TRUE, use.names = TRUE)
  })
  all_names <- unique(unlist(lapply(flattened, names), use.names = FALSE))
  out <- matrix(NA_real_, length(flattened), length(all_names),
                dimnames = list(NULL, all_names))
  for (i in seq_along(flattened)) out[i, names(flattened[[i]])] <- flattened[[i]]
  as.data.frame(out, check.names = FALSE)
}

.MCMfit_dynamic <- function(model, data, weights, compute_se, se_type,
                            optimizers, optim_iters, loss_type,
                            learning_rate, use_bounds, use_skewness,
                            use_kurtosis, device, outofbounds_penalty,
                            monitor_grads, debug, n_starts, seed,
                            stationarity_limit, stationarity_penalty,
                            verbose, moment_weighting, se_correction,
                            weight_ridge, jacobian_method) {
  start_time <- Sys.time()
  if (!isTRUE(use_skewness) || !isTRUE(use_kurtosis)) {
    stop(
      "The dynamic kernel currently requires both `use_skewness = TRUE` and `use_kurtosis = TRUE`.",
      call. = FALSE
    )
  }
  moment_weighting <- .normalize_moment_weighting(moment_weighting)
  se_correction <- .normalize_se_correction(se_correction)
  if (isTRUE(compute_se) && !identical(se_type, "asymptotic")) {
    stop("The dynamic kernel currently supports only `se_type = \"asymptotic\"`.", call. = FALSE)
  }
  if ((isTRUE(compute_se) || !identical(moment_weighting, "identity")) &&
      !identical(loss_type, "mse")) {
    stop(
      "Dynamic WLS weighting and asymptotic SEs require `loss_type = \"mse\"`.",
      call. = FALSE
    )
  }
  if (isTRUE(compute_se) && identical(se_correction, "model_based") &&
      !identical(moment_weighting, "full")) {
    stop(
      "`se_correction = \"model_based\"` requires `moment_weighting = \"full\"`; use robust SEs for identity or diagonal weights.",
      call. = FALSE
    )
  }
  if (length(stationarity_limit) != 1L || !is.finite(stationarity_limit) ||
      stationarity_limit <= 0 || stationarity_limit >= 1) {
    stop("`stationarity_limit` must be a finite number strictly between 0 and 1.", call. = FALSE)
  }
  if (length(stationarity_penalty) != 1L || !is.finite(stationarity_penalty) ||
      stationarity_penalty <= 0) {
    stop("`stationarity_penalty` must be positive.", call. = FALSE)
  }
  if (length(outofbounds_penalty) != 1L || !is.finite(outofbounds_penalty) ||
      outofbounds_penalty < 0) {
    stop("`outofbounds_penalty` must be non-negative.", call. = FALSE)
  }
  if (length(optim_iters) == 1L) optim_iters <- rep(optim_iters, length(optimizers))
  if (length(learning_rate) == 1L) learning_rate <- rep(learning_rate, length(optimizers))
  if (length(optim_iters) != length(optimizers)) {
    stop("`optim_iters` must have length 1 or `length(optimizers)`.", call. = FALSE)
  }
  if (length(learning_rate) != length(optimizers)) {
    stop("`learning_rate` must have length 1 or `length(optimizers)`.", call. = FALSE)
  }

  if (!inherits(data, "mcmdataclass")) {
    data <- MCMdatasummary(
      data, scale_data = model$meta_data$scale_data, weights = weights,
      prep_asymptotic_se = isTRUE(compute_se) ||
        !identical(moment_weighting, "identity"),
      use_skewness = TRUE, use_kurtosis = TRUE
    )
  }
  if (data$meta_data$ncol != model$meta_data$n_phenotypes) {
    stop(
      "Model expected ", model$meta_data$n_phenotypes, " phenotypes but ",
      data$meta_data$ncol, " were found.", call. = FALSE
    )
  }
  if (isTRUE(model$meta_data$weighted) && is.null(weights)) {
    stop("Model was created with weights but no weights were provided.", call. = FALSE)
  }
  weight_spec <- .dynamic_weight_specification(
    data, moment_weighting = moment_weighting,
    weight_ridge = weight_ridge, require_vcov = isTRUE(compute_se)
  )
  model$meta_data$n_obs <- data$meta_data$N
  model$meta_data$kernel <- "dynamic"
  if (is.null(device)) device <- torch_device("cpu")
  dtype <- torch_float64()
  observed <- list(
    M2 = torch_tensor(data$M2, device = device, dtype = dtype),
    M3 = torch_tensor(data$M3, device = device, dtype = dtype),
    M4 = torch_tensor(data$M4, device = device, dtype = dtype)
  )
  masks <- list(
    m2 = .torch_m2m2v_mask(observed$M2, device, dtype),
    m3 = .torch_m3m2v_mask(observed$M3, device, dtype),
    m4 = .torch_m4m2v_mask(observed$M4, device, dtype)
  )
  lossfunc <- .get_lossfunc(loss_type)
  moment_grids <- .dynamic_moment_grids(model$meta_data$n_phenotypes)
  moment_weight <- torch_tensor(
    weight_spec$W, device = device, dtype = dtype
  )
  starts <- .dynamic_random_starts(model, n_starts, seed, stationarity_limit)
  preparation_done <- Sys.time()

  all_runs <- vector("list", length(starts))
  diagnostics <- vector("list", length(starts))
  best <- NULL
  for (start_index in seq_along(starts)) {
    start_model <- .dynamic_model_with_parameters(model, starts[[start_index]])
    if (isTRUE(verbose)) {
      cat(sprintf("Dynamic start %d/%d\n", start_index, length(starts)))
    }
    run <- tryCatch(
      .torch_fit_dynamic(
        start_model, observed, masks, optimizers, optim_iters,
        learning_rate, lossfunc, device, dtype, use_bounds,
        outofbounds_penalty, stationarity_limit, stationarity_penalty,
        moment_weight, moment_grids, loss_type,
        monitor_grads = monitor_grads, verbose = verbose
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "dynamic_fit_error")
    )
    all_runs[[start_index]] <- run
    if (inherits(run, "dynamic_fit_error")) {
      diagnostics[[start_index]] <- data.frame(
        start = start_index, initial_loss = NA_real_, final_loss = NA_real_,
        convergence = 1L, message = run$error, spectral_radius = NA_real_,
        stationary = FALSE, bounds_ok = FALSE, stringsAsFactors = FALSE
      )
      if (isTRUE(verbose)) cat("  failed: ", run$error, "\n", sep = "")
      next
    }
    admissibility <- .dynamic_parameters_admissible(
      model, run$parameters, use_bounds, stationarity_limit
    )
    run$admissibility <- admissibility
    all_runs[[start_index]] <- run
    diagnostics[[start_index]] <- data.frame(
      start = start_index,
      initial_loss = run$initial_loss,
      final_loss = run$loss,
      convergence = if (admissibility$admissible && is.finite(run$loss)) 0L else 1L,
      message = if (admissibility$admissible && is.finite(run$loss)) {
        "Optimizer sequence completed; admissibility checks passed."
      } else {
        "Optimizer sequence completed without an admissible stationary solution."
      },
      spectral_radius = admissibility$spectral_radius,
      stationary = admissibility$stationary,
      bounds_ok = admissibility$bounds_ok,
      stringsAsFactors = FALSE
    )
    if (isTRUE(verbose)) {
      Psi <- .dynamic_tensor_to_matrix(run$predicted$Psi_G)
      cat(sprintf(
        "  final loss=%.8g, convergence=%d, spectral radius=%.6f\n",
        run$loss, diagnostics[[start_index]]$convergence,
        admissibility$spectral_radius
      ))
      cat("  estimated Gaussian covariance:\n")
      print(Psi)
    }
    if (admissibility$admissible && is.finite(run$loss) &&
        (is.null(best) || run$loss < best$loss)) {
      best <- run
      best$start_index <- start_index
    }
  }
  diagnostics_df <- do.call(rbind, diagnostics)
  if (is.null(best)) {
    stop(
      "No optimization start produced a finite, in-bounds stationary solution. ",
      "Inspect starts, bounds, or increase `n_starts`.", call. = FALSE
    )
  }
  optimization_done <- Sys.time()

  fitted_model <- .dynamic_model_with_parameters(model, best$parameters)
  fitted_model$meta_data$kernel <- "dynamic"
  fitted_model$meta_data$stationarity_limit <- stationarity_limit
  predicted <- lapply(
    best$predicted[c("M2", "M3", "M4", "K4", "within_M2", "Psi_G", "L_G")],
    .dynamic_tensor_to_matrix
  )
  B <- .dynamic_tensor_to_matrix(best$predicted$B)
  variables <- fitted_model$meta_data$original_colnames
  dimnames(B) <- list(variables, variables)
  dimnames(predicted$Psi_G) <- list(variables, variables)
  dimnames(predicted$within_M2) <- list(variables, variables)
  tau <- .dynamic_tensor_to_numeric(best$predicted$tau)
  kappa <- .dynamic_tensor_to_numeric(best$predicted$kappa)
  names(tau) <- names(kappa) <- variables
  observed_base <- list(
    M2 = as.matrix(data$M2), M3 = as.matrix(data$M3), M4 = as.matrix(data$M4)
  )
  residuals <- list(
    M2 = observed_base$M2 - predicted$M2,
    M3 = observed_base$M3 - predicted$M3,
    M4 = observed_base$M4 - predicted$M4
  )
  se_started <- Sys.time()
  asymptotic <- if (isTRUE(compute_se)) {
    .dynamic_asymptotic_se(
      fitted_model, data, weight_spec,
      se_correction = se_correction,
      jacobian_method = jacobian_method,
      weight_ridge = weight_ridge
    )
  } else {
    NULL
  }
  se_finished <- Sys.time()
  parameter_se <- if (is.null(asymptotic)) {
    stats::setNames(rep(NA_real_, length(fitted_model$param_names)),
                    fitted_model$param_names)
  } else {
    asymptotic$se
  }
  transition <- .dynamic_transition_table(fitted_model, B)
  transition$se <- unname(parameter_se[transition$label])
  innovation <- data.frame(
    variable = variables,
    variance = rep(1, length(variables)),
    third_cumulant = tau,
    fourth_cumulant = kappa,
    raw_fourth = kappa + 3,
    third_se = unname(parameter_se[paste0("tau_", variables)]),
    fourth_se = unname(parameter_se[paste0("kappa_", variables)]),
    row.names = NULL,
    check.names = FALSE
  )
  gaussian <- if (isTRUE(fitted_model$meta_data$gaussian_residual)) {
    .dynamic_gaussian_table(
      predicted$Psi_G, variables,
      if (is.null(asymptotic)) NULL else asymptotic$gaussian_se
    )
  } else {
    data.frame(
      label = character(), lhs = character(), rhs = character(),
      estimate = numeric(), se = numeric(), stringsAsFactors = FALSE
    )
  }
  dof <- MCMdegreesoffreedom(fitted_model, TRUE, TRUE)
  result_values <- if (isTRUE(compute_se)) {
    rbind(est = best$parameters, se = unname(parameter_se))
  } else {
    matrix(best$parameters, nrow = 1, dimnames = list("est", NULL))
  }
  result_df <- as.data.frame(result_values, check.names = FALSE)
  colnames(result_df) <- fitted_model$param_names
  fitted_model$param_values <- best$parameters
  fitted_model$start_values$set_all(best$parameters)

  gradient_frame <- .dynamic_gradient_history_frame(best$gradient_history)
  convergence <- list(
    code = 0L,
    message = diagnostics_df$message[best$start_index],
    best_start = best$start_index
  )
  finish_time <- Sys.time()
  info <- list(
    version = MCMSEMversion,
    kernel = "dynamic",
    compute_se = isTRUE(compute_se),
    se_type = if (isTRUE(compute_se)) "asymptotic" else NA_character_,
    optim_iters = optim_iters,
    learning_rate = learning_rate,
    use_bounds = use_bounds,
    use_skewness = TRUE,
    use_kurtosis = TRUE,
    jacobian_method = if (isTRUE(compute_se)) jacobian_method else NA_character_,
    debug = debug,
    device = device$type,
    device_se = NA_character_,
    low_memory = FALSE,
    weighted = data$meta_data$weighted,
    loss_type = loss_type,
    moment_weighting = moment_weighting,
    moment_weight_normalization = weight_spec$normalization,
    weight_ridge = weight_ridge,
    se_correction = if (is.null(asymptotic)) NA_character_ else asymptotic$correction,
    jacobian_rank = if (is.null(asymptotic)) NA_integer_ else asymptotic$jacobian_rank,
    jacobian_condition = if (is.null(asymptotic)) NA_real_ else
      asymptotic$jacobian_condition,
    information_condition = if (is.null(asymptotic)) NA_real_ else
      asymptotic$bread_condition,
    optimizers = optimizers,
    n = data$meta_data$N,
    n_starts = length(starts),
    seed = seed,
    stationarity_limit = stationarity_limit,
    stationarity_penalty = stationarity_penalty,
    degrees_of_freedom = dof,
    convergence = convergence,
    start_diagnostics = diagnostics_df,
    all_start_details = all_runs,
    gradient_history = gradient_frame,
    standard_errors_available = isTRUE(compute_se) &&
      all(is.finite(parameter_se)),
    standard_errors_note = if (isTRUE(compute_se)) {
      paste0(
        "Asymptotic ", if (is.null(asymptotic)) NA_character_ else asymptotic$correction,
        " covariance based on the mean-corrected raw central-moment influence functions."
      )
    } else {
      "Set `compute_se = TRUE` to calculate dynamic-kernel asymptotic standard errors."
    }
  )

  mcmresultclass(
    df = result_df,
    loss = best$loss,
    gradients = mcmmultigradienthistoryclass(),
    model = fitted_model$copy(),
    history = list(loss = best$loss_history, all_starts = diagnostics_df),
    runtimes = list(
      Preparation = preparation_done - start_time,
      Optimizer = optimization_done - preparation_done,
      SE = se_finished - se_started,
      Total = finish_time - start_time
    ),
    info = info,
    observed = observed_base,
    predicted = predicted,
    kernel = "dynamic",
    B = B,
    transition_parameters = transition,
    innovation_variances = stats::setNames(rep(1, length(variables)), variables),
    innovation_third = tau,
    innovation_fourth = kappa,
    Psi_G = predicted$Psi_G,
    spectral_radius = best$admissibility$spectral_radius,
    stationary = TRUE,
    convergence = convergence,
    degrees_of_freedom = dof$df,
    n_moments = dof$n_moments,
    start_diagnostics = diagnostics_df,
    residuals = residuals,
    dynamic = list(
      transition = transition,
      innovations = innovation,
      gaussian_covariance = gaussian,
      within_M2 = predicted$within_M2,
      K4 = predicted$K4,
      L_G = predicted$L_G,
      moment_weighting = moment_weighting,
      moment_weight = weight_spec$W,
      moment_vcov = weight_spec$omega,
      asymptotic = asymptotic
    )
  )
}

MCMdynamicloss <- function(model, data, parameters = NULL,
                           loss_type = "mse",
                           moment_weighting = "identity",
                           weight_ridge = 1e-8) {
  if (!identical(.model_kernel(model), "dynamic")) {
    stop("`model` must use `kernel = \"dynamic\"`.", call. = FALSE)
  }
  if (!inherits(data, "mcmdataclass")) {
    moment_weighting <- .normalize_moment_weighting(moment_weighting)
    data <- MCMdatasummary(
      data, scale_data = model$meta_data$scale_data,
      prep_asymptotic_se = !identical(moment_weighting, "identity")
    )
  }
  implied <- .dynamic_implied_moments_base(model, parameters)
  observed <- c(
    .dynamic_unique_moments(data$M2, 2L),
    .dynamic_unique_moments(data$M3, 3L),
    .dynamic_unique_moments(data$M4, 4L)
  )
  predicted <- c(
    .dynamic_unique_moments(implied$M2, 2L),
    .dynamic_unique_moments(implied$M3, 3L),
    .dynamic_unique_moments(implied$M4, 4L)
  )
  residual <- predicted - observed
  moment_weighting <- .normalize_moment_weighting(moment_weighting)
  if (!identical(moment_weighting, "identity")) {
    if (!identical(loss_type, "mse")) {
      stop("Weighted dynamic loss requires `loss_type = \"mse\"`.", call. = FALSE)
    }
    W <- .dynamic_weight_specification(
      data, moment_weighting, weight_ridge, require_vcov = TRUE
    )$W
    return(drop(crossprod(residual, W %*% residual)))
  }
  switch(
    loss_type,
    mse = sum(residual^2),
    l1 = sum(abs(residual)),
    smooth_l1 = sum(ifelse(abs(residual) < 1, residual^2 / 2,
                           abs(residual) - 0.5)),
    stop("`loss_type` must be one of 'mse', 'l1', or 'smooth_l1'.", call. = FALSE)
  )
}
