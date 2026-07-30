#!/usr/bin/env Rscript

# Opt-in recovery and asymptotic-SE calibration for a stationary dynamic model
# with signed-gamma innovations and one additive common-gamma confounder.
#
# Default design:
#   Rscript inst/validation/common_gamma_monte_carlo.R
#
# Fast smoke design:
#   Rscript inst/validation/common_gamma_monte_carlo.R \
#     --repetitions=2 --n=3000 --starts=1 \
#     --rprop-iters=80 --lbfgs-iters=4

parse_arguments <- function(args) {
  settings <- list(
    repetitions = 50L, n = 5000L, starts = 3L, burnin = 220L,
    rprop_iters = 200L, lbfgs_iters = 10L, master_seed = 20260729L,
    moment_weighting = "diagonal", se_correction = "robust",
    output = "validation-output/common-gamma"
  )
  for (argument in args) {
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    name <- gsub("-", "_", pieces[1L])
    if (length(pieces) != 2L || !(name %in% names(settings))) {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
    settings[[name]] <- if (is.numeric(settings[[name]])) {
      as.integer(pieces[2L])
    } else pieces[2L]
  }
  settings
}

settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to validate the MCMSEM working tree.")
}
devtools::load_all(".", quiet = TRUE)

signed_gamma <- function(n, shape, sign = 1) {
  sign * (stats::rgamma(n, shape = shape) - shape) / sqrt(shape)
}

B_true <- matrix(c(0.42, 0.16, -0.11, 0.32), 2L, byrow = TRUE)
innovation_shapes <- c(X = 4, Y = 6)
innovation_signs <- c(X = 1, Y = -1)
gamma_loadings <- c(X = 0.55, Y = -0.35)
gamma_shape <- 3

truth <- c(
  phi_X = B_true[1L, 1L],
  X_lag_to_Y = B_true[2L, 1L],
  Y_lag_to_X = B_true[1L, 2L],
  phi_Y = B_true[2L, 2L],
  shape_X = innovation_shapes[["X"]],
  shape_Y = innovation_shapes[["Y"]],
  tau_X = innovation_signs[["X"]] * 2 / sqrt(innovation_shapes[["X"]]),
  tau_Y = innovation_signs[["Y"]] * 2 / sqrt(innovation_shapes[["Y"]]),
  kappa_X = 6 / innovation_shapes[["X"]],
  kappa_Y = 6 / innovation_shapes[["Y"]],
  loading_Gamma_X = gamma_loadings[["X"]],
  loading_Gamma_Y = gamma_loadings[["Y"]],
  shape_Gamma = gamma_shape
)

residual_truth <- c(
  UVar_X = gamma_loadings[["X"]]^2,
  UCov_XY = prod(gamma_loadings),
  UVar_Y = gamma_loadings[["Y"]]^2
)

simulate_cross_section <- function(n, seed) {
  set.seed(seed)
  state <- matrix(0, n, 2L)
  for (time in seq_len(settings$burnin)) {
    innovations <- cbind(
      signed_gamma(n, innovation_shapes[["X"]], innovation_signs[["X"]]),
      signed_gamma(n, innovation_shapes[["Y"]], innovation_signs[["Y"]])
    )
    state <- state %*% t(B_true) + innovations
  }
  common <- signed_gamma(n, gamma_shape)
  observed <- state + tcrossprod(common, gamma_loadings)
  data <- as.data.frame(observed)
  names(data) <- names(gamma_loadings)
  data
}

make_model <- function(summary_data) {
  model <- MCMmodel(
    summary_data, kernel = "dynamic", residual_family = "common_gamma"
  )
  for (variable in names(innovation_shapes)) {
    shape <- paste0("shape_", variable)
    sign <- paste0("sign_", variable)
    tau <- paste0("tau_", variable)
    kappa <- paste0("kappa_", variable)
    model <- MCMparameter(
      model, shape, "free", start = innovation_shapes[[variable]],
      transform = "positive"
    )
    model <- MCMparameter(
      model, sign, "fixed", value = innovation_signs[[variable]]
    )
    model <- MCMparameter(
      model, tau, "derived",
      expression = stats::as.formula(
        paste0("~ ", sign, " * 2 / sqrt(", shape, ")")
      )
    )
    model <- MCMparameter(
      model, kappa, "derived",
      expression = stats::as.formula(paste0("~ 6 / ", shape))
    )
  }
  for (name in intersect(names(truth), model$param_names)) {
    model <- MCMedit(model, "start", name, truth[[name]])
  }
  model
}

fit_replication <- function(replication) {
  cat(sprintf(
    "[%s] replication %d/%d: simulating N=%d\n",
    format(Sys.time()), replication, settings$repetitions, settings$n
  ))
  data <- simulate_cross_section(
    settings$n, settings$master_seed + replication
  )
  summary_data <- MCMdatasummary(
    data, scale_data = FALSE, prep_asymptotic_se = TRUE,
    use_skewness = TRUE, use_kurtosis = TRUE
  )
  model <- make_model(summary_data)
  true_loss <- MCMdynamicloss(
    model, summary_data, moment_weighting = settings$moment_weighting
  )
  elapsed <- system.time({
    fit <- MCMfit(
      model, summary_data, compute_se = TRUE,
      moment_weighting = settings$moment_weighting,
      se_correction = settings$se_correction,
      optimizers = c("rprop", "lbfgs"),
      optim_iters = c(settings$rprop_iters, settings$lbfgs_iters),
      learning_rate = c(0.01, 0.2), n_starts = settings$starts,
      seed = settings$master_seed + 100000L + replication,
      verbose = FALSE
    )
  })
  estimates <- as.numeric(fit$df["est", names(truth), drop = TRUE])
  standard_errors <- as.numeric(fit$df["se", names(truth), drop = TRUE])
  residual <- fit$dynamic$residual_covariance
  estimates <- c(estimates, stats::setNames(residual$estimate, names(residual_truth)))
  standard_errors <- c(
    standard_errors, stats::setNames(residual$se, names(residual_truth))
  )
  all_truth <- c(truth, residual_truth)
  names(estimates)[seq_along(truth)] <- names(truth)
  names(standard_errors)[seq_along(truth)] <- names(truth)
  cat(sprintf(
    "[%s] replication %d: loss=%.6g, true loss=%.6g, rho=%.4f, rank=%d, elapsed=%.1fs\n",
    format(Sys.time()), replication, fit$loss, true_loss,
    fit$spectral_radius, fit$dynamic$asymptotic$jacobian_rank,
    elapsed[["elapsed"]]
  ))
  data.frame(
    replication = replication, parameter = names(all_truth),
    truth = unname(all_truth), estimate = unname(estimates[names(all_truth)]),
    se = unname(standard_errors[names(all_truth)]),
    covered_95 = abs(estimates[names(all_truth)] - all_truth) <=
      1.96 * standard_errors[names(all_truth)],
    fit_loss = fit$loss, true_loss = true_loss,
    spectral_radius = fit$spectral_radius,
    jacobian_rank = fit$dynamic$asymptotic$jacobian_rank,
    information_condition = fit$info$information_condition,
    admissible_starts = sum(
      fit$start_diagnostics$convergence == 0 &
        fit$start_diagnostics$stationary & fit$start_diagnostics$bounds_ok
    ),
    elapsed_seconds = elapsed[["elapsed"]], check.names = FALSE
  )
}

start_time <- Sys.time()
rows <- vector("list", settings$repetitions)
errors <- character(settings$repetitions)
for (replication in seq_len(settings$repetitions)) {
  rows[[replication]] <- tryCatch(
    fit_replication(replication),
    error = function(error) {
      errors[replication] <<- conditionMessage(error)
      cat(sprintf("replication %d failed: %s\n", replication, errors[replication]))
      NULL
    }
  )
}
results <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
if (is.null(results) || !nrow(results)) stop("Every replication failed.")

calibration <- do.call(rbind, lapply(split(results, results$parameter), function(x) {
  empirical_sd <- stats::sd(x$estimate, na.rm = TRUE)
  data.frame(
    parameter = x$parameter[1L], truth = x$truth[1L],
    mean_estimate = mean(x$estimate, na.rm = TRUE),
    bias = mean(x$estimate - x$truth, na.rm = TRUE),
    empirical_sd = empirical_sd, mean_se = mean(x$se, na.rm = TRUE),
    median_se = stats::median(x$se, na.rm = TRUE),
    se_to_empirical_sd = mean(x$se, na.rm = TRUE) / empirical_sd,
    median_se_to_empirical_sd = stats::median(x$se, na.rm = TRUE) /
      empirical_sd,
    coverage_95 = mean(x$covered_95, na.rm = TRUE),
    successful = sum(is.finite(x$estimate) & is.finite(x$se)),
    row.names = NULL, check.names = FALSE
  )
}))
run_summary <- data.frame(
  requested_repetitions = settings$repetitions,
  successful_fits = length(unique(results$replication)),
  failure_rate = mean(nzchar(errors)), n = settings$n,
  starts = settings$starts,
  total_elapsed_seconds = as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  ), check.names = FALSE
)

dir.create(settings$output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  results, file.path(settings$output, "replication-results.csv"), row.names = FALSE
)
utils::write.csv(
  calibration, file.path(settings$output, "se-calibration.csv"), row.names = FALSE
)
utils::write.csv(
  run_summary, file.path(settings$output, "run-summary.csv"), row.names = FALSE
)
saveRDS(
  list(settings = settings, truth = c(truth, residual_truth),
       results = results, calibration = calibration,
       run_summary = run_summary, errors = errors),
  file.path(settings$output, "validation-results.rds")
)

print(calibration)
print(run_summary)
cat("Validation output: ", normalizePath(settings$output), "\n", sep = "")
