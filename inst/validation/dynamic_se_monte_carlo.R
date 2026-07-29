#!/usr/bin/env Rscript

# Monte Carlo calibration of dynamic-kernel asymptotic standard errors.
#
# Default design (50 replications):
#   Rscript inst/validation/dynamic_se_monte_carlo.R
#
# Fast smoke design:
#   Rscript inst/validation/dynamic_se_monte_carlo.R \
#     --repetitions=2 --n=3000 --starts=1 --rprop-iters=80 --lbfgs-iters=4
#
# The fitted moments are raw central moments. MCMdatasummary() estimates their
# sampling covariance from mean-corrected influence functions. The fit uses
# identity, diagonal, or full WLS weights, and MCMfit() reports robust sandwich
# SEs except that auto mode selects efficient information-matrix SEs for full
# WLS. Psi_G SEs are obtained by the delta method from its Cholesky parameters.

parse_args <- function(args) {
  defaults <- list(
    repetitions = 50L,
    n = 20000L,
    starts = 3L,
    burnin = 250L,
    rprop_iters = 250L,
    lbfgs_iters = 10L,
    moment_weighting = "diagonal",
    se_correction = "auto",
    gaussian_residual = TRUE,
    include_true_start = TRUE,
    output_dir = "validation-output/dynamic-se"
  )
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L) stop("Arguments must have the form --name=value")
    name <- gsub("-", "_", pieces[1])
    if (!(name %in% names(defaults))) stop("Unknown argument: ", pieces[1])
    defaults[[name]] <- if (is.logical(defaults[[name]])) {
      tolower(pieces[2]) %in% c("true", "t", "1", "yes")
    } else if (is.numeric(defaults[[name]])) {
      as.integer(pieces[2])
    } else {
      pieces[2]
    }
  }
  defaults
}

config <- parse_args(commandArgs(trailingOnly = TRUE))

in_source_tree <- file.exists("DESCRIPTION") &&
  identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1]), "MCMSEM")
if (in_source_tree) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install devtools to validate the MCMSEM working tree.")
  }
  devtools::load_all(".", quiet = TRUE)
} else if (!requireNamespace("MCMSEM", quietly = TRUE)) {
  stop("Run from an MCMSEM source tree or install MCMSEM first.")
}

B_true <- matrix(c(0.45, 0.16, -0.12, 0.35), 2, byrow = TRUE)
Psi_true <- matrix(c(0.35, 0.12, 0.12, 0.30), 2)

covariance_to_parameters <- function(Psi) {
  L <- t(chol(Psi))
  c(
    log_sd_G_X = log(L[1, 1]),
    chol_G_Y_X = L[2, 1],
    log_sd_G_Y = log(L[2, 2])
  )
}

truth_raw <- c(
  phi_X = B_true[1, 1],
  X_lag_to_Y = B_true[2, 1],
  Y_lag_to_X = B_true[1, 2],
  phi_Y = B_true[2, 2],
  tau_X = 1,
  tau_Y = -sqrt(8 / 5),
  kappa_X = 1.5,
  kappa_Y = 12 / 5,
  covariance_to_parameters(Psi_true)
)
truth_report_all <- c(
  truth_raw[c("phi_X", "X_lag_to_Y", "Y_lag_to_X", "phi_Y",
              "tau_X", "tau_Y", "kappa_X", "kappa_Y")],
  GVar_X = Psi_true[1, 1],
  GCov_XY = Psi_true[2, 1],
  GVar_Y = Psi_true[2, 2]
)
truth_report <- if (isTRUE(config$gaussian_residual)) {
  truth_report_all
} else {
  truth_report_all[seq_len(8L)]
}

simulate_final_cross_section <- function(n, seed) {
  set.seed(seed)
  state <- matrix(0, n, 2)
  for (time in seq_len(config$burnin + 10L)) {
    innovation <- cbind(
      (stats::rgamma(n, shape = 4) - 4) / 2,
      -(stats::rchisq(n, df = 5) - 5) / sqrt(10)
    )
    state <- state %*% t(B_true) + innovation
  }
  observed <- if (isTRUE(config$gaussian_residual)) {
    state + MASS::mvrnorm(n, mu = c(0, 0), Sigma = Psi_true)
  } else {
    state
  }
  data <- as.data.frame(observed)
  names(data) <- c("X", "Y")
  data
}

fit_replication <- function(replication) {
  cat(sprintf("[%s] replication %d/%d\n", format(Sys.time()),
              replication, config$repetitions))
  data <- simulate_final_cross_section(config$n, 810000L + replication)
  summary_data <- MCMSEM::MCMdatasummary(
    data, scale_data = FALSE, prep_asymptotic_se = TRUE,
    use_skewness = TRUE, use_kurtosis = TRUE
  )
  model <- MCMSEM::MCMmodel(
    summary_data, kernel = "dynamic",
    gaussian_residual = config$gaussian_residual
  )
  if (isTRUE(config$include_true_start)) {
    for (name in intersect(names(truth_raw), model$param_names)) {
      model <- MCMSEM::MCMedit(model, "start", name, truth_raw[[name]])
    }
  }
  elapsed <- system.time({
    fit <- MCMSEM::MCMfit(
      model, summary_data,
      compute_se = TRUE, se_type = "asymptotic",
      moment_weighting = config$moment_weighting,
      se_correction = config$se_correction,
      optimizers = c("rprop", "lbfgs"),
      optim_iters = c(config$rprop_iters, config$lbfgs_iters),
      learning_rate = c(0.01, 0.2),
      n_starts = config$starts,
      seed = 850000L + replication,
      verbose = FALSE
    )
  })

  estimates <- c(
    phi_X = fit$df["est", "phi_X"],
    X_lag_to_Y = fit$df["est", "X_lag_to_Y"],
    Y_lag_to_X = fit$df["est", "Y_lag_to_X"],
    phi_Y = fit$df["est", "phi_Y"],
    tau_X = fit$df["est", "tau_X"],
    tau_Y = fit$df["est", "tau_Y"],
    kappa_X = fit$df["est", "kappa_X"],
    kappa_Y = fit$df["est", "kappa_Y"]
  )
  standard_errors <- c(
    phi_X = fit$df["se", "phi_X"],
    X_lag_to_Y = fit$df["se", "X_lag_to_Y"],
    Y_lag_to_X = fit$df["se", "Y_lag_to_X"],
    phi_Y = fit$df["se", "phi_Y"],
    tau_X = fit$df["se", "tau_X"],
    tau_Y = fit$df["se", "tau_Y"],
    kappa_X = fit$df["se", "kappa_X"],
    kappa_Y = fit$df["se", "kappa_Y"]
  )
  if (isTRUE(config$gaussian_residual)) {
    gaussian <- fit$dynamic$gaussian_covariance
    estimates <- c(
      estimates, GVar_X = gaussian$estimate[1],
      GCov_XY = gaussian$estimate[2], GVar_Y = gaussian$estimate[3]
    )
    standard_errors <- c(
      standard_errors, GVar_X = gaussian$se[1],
      GCov_XY = gaussian$se[2], GVar_Y = gaussian$se[3]
    )
  }
  long <- data.frame(
    replication = replication,
    parameter = names(truth_report),
    truth = unname(truth_report),
    estimate = unname(estimates[names(truth_report)]),
    se = unname(standard_errors[names(truth_report)]),
    covered_95 = abs(estimates[names(truth_report)] - truth_report) <=
      1.96 * standard_errors[names(truth_report)],
    fit_loss = fit$loss,
    spectral_radius = fit$spectral_radius,
    convergence = fit$convergence$code,
    jacobian_rank = fit$dynamic$asymptotic$jacobian_rank,
    elapsed_seconds = elapsed[["elapsed"]],
    check.names = FALSE
  )
  cat(sprintf(
    "  loss=%.6g, rho=%.4f, rank=%d, median(SE)=%.5g, elapsed=%.1fs\n",
    fit$loss, fit$spectral_radius, fit$dynamic$asymptotic$jacobian_rank,
    stats::median(standard_errors), elapsed[["elapsed"]]
  ))
  long
}

start_time <- Sys.time()
rows <- vector("list", config$repetitions)
errors <- character(config$repetitions)
for (replication in seq_len(config$repetitions)) {
  rows[[replication]] <- tryCatch(
    fit_replication(replication),
    error = function(e) {
      errors[replication] <<- conditionMessage(e)
      cat(sprintf("  FAILED: %s\n", errors[replication]))
      NULL
    }
  )
}
results <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
if (is.null(results) || nrow(results) == 0L) stop("Every replication failed.")

calibration <- do.call(rbind, lapply(split(results, results$parameter), function(x) {
  data.frame(
    parameter = x$parameter[1],
    truth = x$truth[1],
    mean_estimate = mean(x$estimate, na.rm = TRUE),
    bias = mean(x$estimate - x$truth, na.rm = TRUE),
    empirical_sd = stats::sd(x$estimate, na.rm = TRUE),
    mean_se = mean(x$se, na.rm = TRUE),
    median_se = stats::median(x$se, na.rm = TRUE),
    se_to_empirical_sd = mean(x$se, na.rm = TRUE) /
      stats::sd(x$estimate, na.rm = TRUE),
    coverage_95 = mean(x$covered_95, na.rm = TRUE),
    successful = sum(is.finite(x$estimate) & is.finite(x$se)),
    row.names = NULL,
    check.names = FALSE
  )
}))
run_summary <- data.frame(
  requested_repetitions = config$repetitions,
  successful_fits = length(unique(results$replication)),
  failure_rate = mean(nzchar(errors)),
  mean_elapsed_seconds = mean(results$elapsed_seconds, na.rm = TRUE),
  total_elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  moment_weighting = config$moment_weighting,
  se_correction = config$se_correction,
  n = config$n,
  starts = config$starts,
  gaussian_residual = config$gaussian_residual,
  check.names = FALSE
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, file.path(config$output_dir, "replication-results.csv"), row.names = FALSE)
utils::write.csv(calibration, file.path(config$output_dir, "se-calibration.csv"), row.names = FALSE)
utils::write.csv(run_summary, file.path(config$output_dir, "run-summary.csv"), row.names = FALSE)
saveRDS(
  list(config = config, truth = truth_report, results = results,
       calibration = calibration, run_summary = run_summary, errors = errors),
  file.path(config$output_dir, "validation-results.rds")
)

print(calibration)
print(run_summary)
cat("Validation output: ", normalizePath(config$output_dir), "\n", sep = "")
