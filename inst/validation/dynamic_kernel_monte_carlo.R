#!/usr/bin/env Rscript

# Manual validation for Stationary Dynamic MCMSEM.
#
# Full design requested for methodological validation:
#   Rscript inst/validation/dynamic_kernel_monte_carlo.R
#
# Fast smoke run:
#   Rscript inst/validation/dynamic_kernel_monte_carlo.R \
#     --repetitions=1 --n=5000 --starts=2 --rprop-iters=150 --lbfgs-iters=8
#
# Results are written only when this script is run explicitly. Generated files
# belong in validation-output/, which is excluded from package builds and Git.

parse_args <- function(args) {
  defaults <- list(
    repetitions = 20L,
    n = 200000L,
    starts = 5L,
    burnin = 250L,
    rprop_iters = 300L,
    lbfgs_iters = 12L,
    output_dir = "validation-output/dynamic-kernel",
    include_true_start = TRUE
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
innovation_truth <- c(
  tau_X = 1,
  tau_Y = -sqrt(8 / 5),
  kappa_X = 1.5,
  kappa_Y = 12 / 5
)

covariance_to_parameters <- function(Psi) {
  L <- t(chol(Psi))
  c(
    log_sd_G_X = log(L[1, 1]),
    chol_G_Y_X = L[2, 1],
    log_sd_G_Y = log(L[2, 2])
  )
}

truth <- c(
  phi_X = B_true[1, 1],
  X_lag_to_Y = B_true[2, 1],
  Y_lag_to_X = B_true[1, 2],
  phi_Y = B_true[2, 2],
  innovation_truth,
  covariance_to_parameters(Psi_true)
)

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
  gaussian <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Psi_true)
  out <- as.data.frame(state + gaussian)
  names(out) <- c("X", "Y")
  out
}

fit_replication <- function(replication) {
  cat(sprintf(
    "[%s] replication %d/%d: simulating N=%d\n",
    format(Sys.time()), replication, config$repetitions, config$n
  ))
  data <- simulate_final_cross_section(config$n, 900000L + replication)
  summary_data <- MCMSEM::MCMdatasummary(
    data, scale_data = FALSE, prep_asymptotic_se = FALSE
  )
  model <- MCMSEM::MCMmodel(summary_data, kernel = "dynamic")
  true_loss <- MCMSEM::MCMdynamicloss(model, summary_data, truth)
  if (isTRUE(config$include_true_start)) {
    for (name in names(truth)) {
      model <- MCMSEM::MCMedit(model, "start", name, truth[[name]])
    }
  }
  elapsed <- system.time({
    fit <- MCMSEM::MCMfit(
      model, summary_data,
      compute_se = FALSE,
      optimizers = c("rprop", "lbfgs"),
      optim_iters = c(config$rprop_iters, config$lbfgs_iters),
      learning_rate = c(0.01, 0.2),
      n_starts = config$starts,
      seed = 950000L + replication,
      verbose = FALSE
    )
  })
  estimates <- c(
    phi_X = fit$B[1, 1],
    X_lag_to_Y = fit$B[2, 1],
    Y_lag_to_X = fit$B[1, 2],
    phi_Y = fit$B[2, 2],
    tau_X = fit$innovation_third[["X"]],
    tau_Y = fit$innovation_third[["Y"]],
    kappa_X = fit$innovation_fourth[["X"]],
    kappa_Y = fit$innovation_fourth[["Y"]],
    GVar_X = fit$Psi_G[1, 1],
    GCov_XY = fit$Psi_G[1, 2],
    GVar_Y = fit$Psi_G[2, 2]
  )
  cat(sprintf(
    "[%s] replication %d: fit loss=%.7g, true loss=%.7g, rho=%.5f, elapsed=%.1fs\n",
    format(Sys.time()), replication, fit$loss, true_loss,
    fit$spectral_radius, elapsed[["elapsed"]]
  ))
  data.frame(
    replication = replication,
    as.list(estimates),
    fit_loss = fit$loss,
    true_loss = true_loss,
    convergence = fit$convergence$code,
    stationary = fit$stationary,
    best_start = fit$convergence$best_start,
    elapsed_seconds = elapsed[["elapsed"]],
    check.names = FALSE
  )
}

start_time <- Sys.time()
rows <- vector("list", config$repetitions)
errors <- character(config$repetitions)
for (replication in seq_len(config$repetitions)) {
  rows[[replication]] <- tryCatch(
    fit_replication(replication),
    error = function(e) {
      errors[replication] <<- conditionMessage(e)
      cat(sprintf("replication %d failed: %s\n", replication, errors[replication]))
      data.frame(replication = replication, fit_loss = NA_real_,
                 true_loss = NA_real_, convergence = 1L,
                 stationary = FALSE, best_start = NA_integer_,
                 elapsed_seconds = NA_real_)
    }
  )
}
estimates <- do.call(rbind, rows)

truth_report <- c(
  phi_X = truth[["phi_X"]], X_lag_to_Y = truth[["X_lag_to_Y"]],
  Y_lag_to_X = truth[["Y_lag_to_X"]], phi_Y = truth[["phi_Y"]],
  tau_X = truth[["tau_X"]], tau_Y = truth[["tau_Y"]],
  kappa_X = truth[["kappa_X"]], kappa_Y = truth[["kappa_Y"]],
  GVar_X = Psi_true[1, 1], GCov_XY = Psi_true[1, 2],
  GVar_Y = Psi_true[2, 2]
)
parameter_names <- intersect(names(truth_report), names(estimates))
summary <- do.call(rbind, lapply(parameter_names, function(name) {
  estimate <- estimates[[name]]
  target <- truth_report[[name]]
  data.frame(
    parameter = name,
    truth = target,
    mean = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - target, na.rm = TRUE),
    sd = stats::sd(estimate, na.rm = TRUE),
    median = stats::median(estimate, na.rm = TRUE),
    rmse = sqrt(mean((estimate - target)^2, na.rm = TRUE)),
    successful = sum(is.finite(estimate)),
    check.names = FALSE
  )
}))
run_summary <- data.frame(
  requested_repetitions = config$repetitions,
  successful_fits = sum(is.finite(estimates$fit_loss)),
  failure_rate = mean(!is.finite(estimates$fit_loss)),
  convergence_rate = mean(estimates$convergence == 0, na.rm = TRUE),
  stationary_rate = mean(estimates$stationary, na.rm = TRUE),
  mean_fit_loss = mean(estimates$fit_loss, na.rm = TRUE),
  mean_true_loss = mean(estimates$true_loss, na.rm = TRUE),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(estimates, file.path(config$output_dir, "estimates.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(config$output_dir, "parameter-summary.csv"), row.names = FALSE)
utils::write.csv(run_summary, file.path(config$output_dir, "run-summary.csv"), row.names = FALSE)
saveRDS(
  list(config = config, truth = truth_report, estimates = estimates,
       parameter_summary = summary, run_summary = run_summary, errors = errors),
  file.path(config$output_dir, "validation-results.rds")
)

print(summary)
print(run_summary)
cat("Validation output: ", normalizePath(config$output_dir), "\n", sep = "")
