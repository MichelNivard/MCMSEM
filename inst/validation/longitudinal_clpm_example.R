#!/usr/bin/env Rscript

# Reproduce the longitudinal/largest-wave comparisons reported in README.md,
# wiki/2.3 Choosing a kernel.md, and longitudinal_clpm_example.md.
#
# Run from the package source tree:
#   Rscript inst/validation/longitudinal_clpm_example.R

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Install the suggested package `lavaan` to run this example.")
}

in_source_tree <- file.exists("DESCRIPTION") &&
  identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1]), "MCMSEM")
if (in_source_tree) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Install `devtools` to validate the MCMSEM working tree.")
  }
  devtools::load_all(".", quiet = TRUE)
} else if (!requireNamespace("MCMSEM", quietly = TRUE)) {
  stop("Run from an MCMSEM source tree or install MCMSEM first.")
}
arguments <- commandArgs(trailingOnly = TRUE)
run_unrestricted_gamma <- "--unrestricted-gamma" %in% arguments
run_simulation_only <- "--simulation-only" %in% arguments

extract_dynamic_paths <- function(fit) {
  paths <- fit$transition_parameters[
    , c("label", "lagged", "current", "estimate", "se")
  ]
  paths$pvalue <- 2 * stats::pnorm(
    abs(paths$estimate / paths$se), lower.tail = FALSE
  )
  paths$ci.lower <- paths$estimate - 1.96 * paths$se
  paths$ci.upper <- paths$estimate + 1.96 * paths$se
  paths
}

# First validate the comparison under known, exactly matched assumptions.
set.seed(20260730)
simulation_n <- 20000L
simulation_truth <- matrix(
  c(0.55, 0.16, -0.12, 0.45), 2, 2, byrow = TRUE,
  dimnames = list(current = c("X", "Y"), lagged = c("X", "Y"))
)
simulation_state <- matrix(0, simulation_n, 2)
simulation_panel <- array(NA_real_, dim = c(simulation_n, 2, 4))
simulation_innovation <- function(n) {
  cbind(
    stats::rexp(n) - 1,
    (stats::rchisq(n, df = 5) - 5) / sqrt(10)
  )
}
for (tt in seq_len(204L)) {
  simulation_state <- simulation_state %*% t(simulation_truth) +
    simulation_innovation(simulation_n)
  if (tt > 200L) {
    simulation_panel[, , tt - 200L] <- simulation_state
  }
}
simulation_wide <- data.frame(
  X1 = simulation_panel[, 1, 1], Y1 = simulation_panel[, 2, 1],
  X2 = simulation_panel[, 1, 2], Y2 = simulation_panel[, 2, 2],
  X3 = simulation_panel[, 1, 3], Y3 = simulation_panel[, 2, 3],
  X4 = simulation_panel[, 1, 4], Y4 = simulation_panel[, 2, 4]
)
simulation_clpm_syntax <- "
  X2 ~ x_ar*X1 + y_to_x*Y1
  X3 ~ x_ar*X2 + y_to_x*Y2
  X4 ~ x_ar*X3 + y_to_x*Y3
  Y2 ~ y_ar*Y1 + x_to_y*X1
  Y3 ~ y_ar*Y2 + x_to_y*X2
  Y4 ~ y_ar*Y3 + x_to_y*X3
  X1 ~~ Y1
  X2 ~~ Y2
  X3 ~~ Y3
  X4 ~~ Y4
"
simulation_clpm_fit <- lavaan::sem(
  simulation_clpm_syntax, data = simulation_wide,
  estimator = "MLR", meanstructure = TRUE
)
simulation_labels <- c("x_ar", "y_to_x", "x_to_y", "y_ar")
simulation_clpm_estimates <- lavaan::parameterEstimates(
  simulation_clpm_fit, ci = TRUE
)
simulation_clpm_paths <- simulation_clpm_estimates[
  match(simulation_labels, simulation_clpm_estimates$label),
  c("label", "est", "se", "pvalue", "ci.lower", "ci.upper")
]
simulation_final <- simulation_wide[c("X4", "Y4")]
names(simulation_final) <- c("X", "Y")
simulation_data <- MCMdatasummary(
  simulation_final,
  scale_data = FALSE,
  prep_asymptotic_se = TRUE,
  use_skewness = TRUE,
  use_kurtosis = TRUE
)
simulation_model <- MCMmodel(
  simulation_data,
  n_latent = 0,
  kernel = "dynamic",
  gaussian_residual = FALSE
)
simulation_model <- MCMedit(simulation_model, "B", c(1, 1), "x_ar")
simulation_model <- MCMedit(simulation_model, "B", c(1, 2), "y_to_x")
simulation_model <- MCMedit(simulation_model, "B", c(2, 1), "x_to_y")
simulation_model <- MCMedit(simulation_model, "B", c(2, 2), "y_ar")
simulation_dynamic_fit <- MCMfit(
  simulation_model,
  simulation_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(750, 40),
  learning_rate = c(0.01, 0.2),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 20,
  seed = 20260730,
  verbose = FALSE
)
simulation_dynamic_paths <- extract_dynamic_paths(simulation_dynamic_fit)
simulation_metrics <- data.frame(
  n = simulation_n,
  clpm_cfi_robust = lavaan::fitMeasures(simulation_clpm_fit, "cfi.robust"),
  clpm_tli_robust = lavaan::fitMeasures(simulation_clpm_fit, "tli.robust"),
  clpm_rmsea_robust = lavaan::fitMeasures(simulation_clpm_fit, "rmsea.robust"),
  clpm_srmr = lavaan::fitMeasures(simulation_clpm_fit, "srmr"),
  dynamic_loss = simulation_dynamic_fit$loss,
  dynamic_spectral_radius = simulation_dynamic_fit$spectral_radius,
  dynamic_nominal_df = simulation_dynamic_fit$degrees_of_freedom,
  dynamic_jacobian_rank = simulation_dynamic_fit$info$jacobian_rank,
  dynamic_information_condition =
    simulation_dynamic_fit$info$information_condition,
  dynamic_admissible_starts = sum(
    simulation_dynamic_fit$start_diagnostics$convergence == 0
  )
)

# Simulate a second stationary process with centered gamma innovations and add
# a stable bivariate Gaussian random intercept. The innovation shape constraints
# and random-intercept distribution are exactly matched in the fitted models.
simulation_gaussian_n <- 100000L
simulation_shape_x <- 1
simulation_shape_y <- 2.5
simulation_gaussian_truth <- matrix(
  c(0.50, 0.20, 0.20, 0.35), 2, 2, byrow = TRUE,
  dimnames = list(c("X", "Y"), c("X", "Y"))
)
set.seed(20260732)
simulation_ri_state <- matrix(0, simulation_gaussian_n, 2L)
simulation_ri_panel <- array(
  NA_real_, dim = c(simulation_gaussian_n, 2L, 4L)
)
simulation_ri_innovation <- function(n) {
  cbind(
    (stats::rgamma(n, shape = simulation_shape_x) -
       simulation_shape_x) / sqrt(simulation_shape_x),
    (stats::rgamma(n, shape = simulation_shape_y) -
       simulation_shape_y) / sqrt(simulation_shape_y)
  )
}
for (tt in seq_len(204L)) {
  simulation_ri_state <-
    simulation_ri_state %*% t(simulation_truth) +
    simulation_ri_innovation(simulation_gaussian_n)
  if (tt > 200L) {
    simulation_ri_panel[, , tt - 200L] <- simulation_ri_state
  }
}
simulation_random_intercept <- matrix(
  stats::rnorm(simulation_gaussian_n * 2L), simulation_gaussian_n, 2L
) %*% chol(simulation_gaussian_truth)
for (wave in seq_len(4L)) {
  simulation_ri_panel[, , wave] <-
    simulation_ri_panel[, , wave] + simulation_random_intercept
}
simulation_ri_wide <- data.frame(
  X1 = simulation_ri_panel[, 1, 1], Y1 = simulation_ri_panel[, 2, 1],
  X2 = simulation_ri_panel[, 1, 2], Y2 = simulation_ri_panel[, 2, 2],
  X3 = simulation_ri_panel[, 1, 3], Y3 = simulation_ri_panel[, 2, 3],
  X4 = simulation_ri_panel[, 1, 4], Y4 = simulation_ri_panel[, 2, 4]
)

simulation_confounded_clpm_fit <- lavaan::sem(
  simulation_clpm_syntax, data = simulation_ri_wide,
  estimator = "MLR", meanstructure = TRUE
)
simulation_confounded_clpm_estimates <- lavaan::parameterEstimates(
  simulation_confounded_clpm_fit, ci = TRUE
)
simulation_confounded_clpm_paths <- simulation_confounded_clpm_estimates[
  match(simulation_labels, simulation_confounded_clpm_estimates$label),
  c("label", "est", "se", "pvalue", "ci.lower", "ci.upper")
]

simulation_riclpm_syntax <- "
  RI_X =~ 1*X1 + 1*X2 + 1*X3 + 1*X4
  RI_Y =~ 1*Y1 + 1*Y2 + 1*Y3 + 1*Y4
  wX1 =~ 1*X1
  wX2 =~ 1*X2
  wX3 =~ 1*X3
  wX4 =~ 1*X4
  wY1 =~ 1*Y1
  wY2 =~ 1*Y2
  wY3 =~ 1*Y3
  wY4 =~ 1*Y4
  X1 ~~ 0*X1
  X2 ~~ 0*X2
  X3 ~~ 0*X3
  X4 ~~ 0*X4
  Y1 ~~ 0*Y1
  Y2 ~~ 0*Y2
  Y3 ~~ 0*Y3
  Y4 ~~ 0*Y4
  wX2 ~ x_ar*wX1 + y_to_x*wY1
  wX3 ~ x_ar*wX2 + y_to_x*wY2
  wX4 ~ x_ar*wX3 + y_to_x*wY3
  wY2 ~ y_ar*wY1 + x_to_y*wX1
  wY3 ~ y_ar*wY2 + x_to_y*wX2
  wY4 ~ y_ar*wY3 + x_to_y*wX3
  RI_X ~~ RI_Y
  wX1 ~~ wY1
  wX2 ~~ wY2
  wX3 ~~ wY3
  wX4 ~~ wY4
  RI_X ~~ 0*wX1 + 0*wY1
  RI_Y ~~ 0*wX1 + 0*wY1
"
simulation_riclpm_fit <- lavaan::sem(
  simulation_riclpm_syntax, data = simulation_ri_wide,
  estimator = "MLR", meanstructure = TRUE, fixed.x = FALSE
)
if (!isTRUE(lavaan::lavInspect(simulation_riclpm_fit, "post.check"))) {
  stop("The simulated RI-CLPM solution failed lavaan's post-estimation check.")
}
simulation_riclpm_estimates <- lavaan::parameterEstimates(
  simulation_riclpm_fit, ci = TRUE
)
simulation_riclpm_paths <- simulation_riclpm_estimates[
  match(simulation_labels, simulation_riclpm_estimates$label),
  c("label", "est", "se", "pvalue", "ci.lower", "ci.upper")
]

simulation_ri_final <- simulation_ri_wide[c("X4", "Y4")]
names(simulation_ri_final) <- c("X", "Y")
simulation_ri_data <- MCMdatasummary(
  simulation_ri_final,
  scale_data = FALSE,
  prep_asymptotic_se = TRUE,
  use_skewness = TRUE,
  use_kurtosis = TRUE
)
simulation_gaussian_model <- MCMmodel(
  simulation_ri_data,
  n_latent = 0,
  kernel = "dynamic",
  residual_family = "gaussian"
)
simulation_gaussian_model <- MCMedit(
  simulation_gaussian_model, "B", c(1, 1), "x_ar"
)
simulation_gaussian_model <- MCMedit(
  simulation_gaussian_model, "B", c(1, 2), "y_to_x"
)
simulation_gaussian_model <- MCMedit(
  simulation_gaussian_model, "B", c(2, 1), "x_to_y"
)
simulation_gaussian_model <- MCMedit(
  simulation_gaussian_model, "B", c(2, 2), "y_ar"
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "shape_X", "free",
  start = simulation_shape_x, transform = "positive"
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "shape_Y", "free",
  start = simulation_shape_y, transform = "positive"
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "sign_X", "fixed", value = 1
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "sign_Y", "fixed", value = 1
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "tau_X", "derived",
  expression = ~ sign_X * 2 / sqrt(shape_X)
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "kappa_X", "derived",
  expression = ~ 6 / shape_X
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "tau_Y", "derived",
  expression = ~ sign_Y * 2 / sqrt(shape_Y)
)
simulation_gaussian_model <- MCMparameter(
  simulation_gaussian_model, "kappa_Y", "derived",
  expression = ~ 6 / shape_Y
)
simulation_gaussian_cholesky <-
  MCMSEM:::.covariance_to_cholesky_parameters(simulation_gaussian_truth)
simulation_gaussian_starts <- c(
  x_ar = simulation_truth[1, 1],
  y_to_x = simulation_truth[1, 2],
  x_to_y = simulation_truth[2, 1],
  y_ar = simulation_truth[2, 2],
  shape_X = simulation_shape_x,
  shape_Y = simulation_shape_y,
  log_sd_G_X = simulation_gaussian_cholesky[1, 1],
  chol_G_Y_X = simulation_gaussian_cholesky[2, 1],
  log_sd_G_Y = simulation_gaussian_cholesky[2, 2]
)
for (name in names(simulation_gaussian_starts)) {
  simulation_gaussian_model <- MCMedit(
    simulation_gaussian_model, "start", name,
    simulation_gaussian_starts[[name]]
  )
}
simulation_gaussian_fit <- MCMfit(
  simulation_gaussian_model,
  simulation_ri_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(700, 60),
  learning_rate = c(0.01, 0.05),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 10,
  seed = 20260732,
  verbose = FALSE
)
simulation_gaussian_paths <- extract_dynamic_paths(simulation_gaussian_fit)
simulation_gaussian_metrics <- data.frame(
  n = simulation_gaussian_n,
  clpm_cfi_robust = lavaan::fitMeasures(
    simulation_confounded_clpm_fit, "cfi.robust"
  ),
  clpm_rmsea_robust = lavaan::fitMeasures(
    simulation_confounded_clpm_fit, "rmsea.robust"
  ),
  clpm_srmr = lavaan::fitMeasures(
    simulation_confounded_clpm_fit, "srmr"
  ),
  riclpm_cfi_robust = lavaan::fitMeasures(
    simulation_riclpm_fit, "cfi.robust"
  ),
  riclpm_rmsea_robust = lavaan::fitMeasures(
    simulation_riclpm_fit, "rmsea.robust"
  ),
  riclpm_srmr = lavaan::fitMeasures(simulation_riclpm_fit, "srmr"),
  dynamic_loss = simulation_gaussian_fit$loss,
  dynamic_spectral_radius = simulation_gaussian_fit$spectral_radius,
  dynamic_nominal_df = simulation_gaussian_fit$degrees_of_freedom,
  dynamic_jacobian_rank = simulation_gaussian_fit$info$jacobian_rank,
  dynamic_information_condition =
    simulation_gaussian_fit$info$information_condition,
  dynamic_admissible_starts = sum(
    simulation_gaussian_fit$start_diagnostics$convergence == 0
  )
)

if (run_simulation_only) {
  cat("Simulation truth\n")
  print(simulation_truth)
  cat("\nNo-confounder simulation paths\n")
  print(simulation_clpm_paths, row.names = FALSE)
  print(simulation_dynamic_paths, row.names = FALSE)
  cat("\nGaussian random-intercept truth\n")
  print(simulation_gaussian_truth)
  cat("\nGaussian-confounder simulation metrics\n")
  print(simulation_gaussian_metrics, row.names = FALSE)
  cat("\nConfounded CLPM paths\n")
  print(simulation_confounded_clpm_paths, row.names = FALSE)
  cat("\nRI-CLPM paths\n")
  print(simulation_riclpm_paths, row.names = FALSE)
  cat("\nGaussian-residual MCMSEM paths\n")
  print(simulation_gaussian_paths, row.names = FALSE)
  cat("\nEstimated Gaussian residual covariance\n")
  print(simulation_gaussian_fit$Psi_G)
  quit(save = "no", status = 0L)
}

# Then run the comparison on the small, derived SIPP matrix bundled with the
# package. The original Census wave files are not included in the repository.
data_file <- if (in_source_tree) {
  file.path("inst", "extdata", "sipp_2014_panel.csv.gz")
} else {
  system.file("extdata", "sipp_2014_panel.csv.gz", package = "MCMSEM")
}
if (!nzchar(data_file) || !file.exists(data_file)) {
  stop("Cannot find the bundled SIPP analysis matrix.")
}
wide <- utils::read.csv(data_file, na.strings = c("", "NA"))
expected_names <- as.vector(rbind(
  paste0("Earnings", 1:4), paste0("Hours", 1:4)
))
stopifnot(identical(names(wide), expected_names), nrow(wide) == 24505L)

clpm_syntax <- "
  Earnings2 ~ earnings_ar*Earnings1 + hours_to_earnings*Hours1
  Earnings3 ~ earnings_ar*Earnings2 + hours_to_earnings*Hours2
  Earnings4 ~ earnings_ar*Earnings3 + hours_to_earnings*Hours3
  Hours2 ~ hours_ar*Hours1 + earnings_to_hours*Earnings1
  Hours3 ~ hours_ar*Hours2 + earnings_to_hours*Earnings2
  Hours4 ~ hours_ar*Hours3 + earnings_to_hours*Earnings3
  Earnings1 ~~ Hours1
  Earnings2 ~~ Hours2
  Earnings3 ~~ Hours3
  Earnings4 ~~ Hours4
"
clpm_fit <- lavaan::sem(
  clpm_syntax, data = wide,
  estimator = "MLR", missing = "fiml", meanstructure = TRUE
)

riclpm_syntax <- "
  RI_Earnings =~ 1*Earnings1 + 1*Earnings2 + 1*Earnings3 + 1*Earnings4
  RI_Hours =~ 1*Hours1 + 1*Hours2 + 1*Hours3 + 1*Hours4
  wE1 =~ 1*Earnings1
  wE2 =~ 1*Earnings2
  wE3 =~ 1*Earnings3
  wE4 =~ 1*Earnings4
  wH1 =~ 1*Hours1
  wH2 =~ 1*Hours2
  wH3 =~ 1*Hours3
  wH4 =~ 1*Hours4
  Earnings1 ~~ 0*Earnings1
  Earnings2 ~~ 0*Earnings2
  Earnings3 ~~ 0*Earnings3
  Earnings4 ~~ 0*Earnings4
  Hours1 ~~ 0*Hours1
  Hours2 ~~ 0*Hours2
  Hours3 ~~ 0*Hours3
  Hours4 ~~ 0*Hours4
  wE2 ~ earnings_ar*wE1 + hours_to_earnings*wH1
  wE3 ~ earnings_ar*wE2 + hours_to_earnings*wH2
  wE4 ~ earnings_ar*wE3 + hours_to_earnings*wH3
  wH2 ~ hours_ar*wH1 + earnings_to_hours*wE1
  wH3 ~ hours_ar*wH2 + earnings_to_hours*wE2
  wH4 ~ hours_ar*wH3 + earnings_to_hours*wE3
  RI_Earnings ~~ RI_Hours
  wE1 ~~ wH1
  wE2 ~~ wH2
  wE3 ~~ wH3
  wE4 ~~ wH4
  RI_Earnings ~~ 0*wE1 + 0*wH1
  RI_Hours ~~ 0*wE1 + 0*wH1
"
riclpm_fit <- lavaan::sem(
  riclpm_syntax, data = wide,
  estimator = "MLR", missing = "fiml", meanstructure = TRUE,
  fixed.x = FALSE
)
if (!isTRUE(lavaan::lavInspect(riclpm_fit, "post.check"))) {
  stop("The RI-CLPM solution failed lavaan's post-estimation check.")
}

path_labels <- c(
  "earnings_ar", "hours_to_earnings", "earnings_to_hours", "hours_ar"
)
extract_lavaan_paths <- function(fit) {
  estimates <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  paths <- estimates[
    !duplicated(estimates$label) & estimates$label %in% path_labels,
    c("label", "est", "se", "pvalue", "std.all", "ci.lower", "ci.upper")
  ]
  paths[match(path_labels, paths$label), ]
}
clpm_paths <- extract_lavaan_paths(clpm_fit)
riclpm_paths <- extract_lavaan_paths(riclpm_fit)

wave_diagnostics <- do.call(rbind, lapply(seq_len(4L), function(wave) {
  values <- stats::na.omit(wide[paste0(c("Earnings", "Hours"), wave)])
  names(values) <- c("Earnings", "Hours")
  data.frame(
    wave = wave, year = 2012L + wave, complete_n = nrow(values),
    earnings_mean = mean(values$Earnings),
    hours_mean = mean(values$Hours),
    earnings_variance = stats::var(values$Earnings),
    hours_variance = stats::var(values$Hours),
    covariance = stats::cov(values$Earnings, values$Hours)
  )
}))
mcm_wave <- wave_diagnostics$wave[which.max(wave_diagnostics$complete_n)]
mcm_values <- stats::na.omit(
  wide[paste0(c("Earnings", "Hours"), mcm_wave)]
)
names(mcm_values) <- c("Earnings", "Hours")
dynamic_data <- MCMdatasummary(
  mcm_values, scale_data = FALSE, prep_asymptotic_se = TRUE,
  use_skewness = TRUE, use_kurtosis = TRUE
)

make_dynamic_model <- function(residual_family) {
  model <- MCMmodel(
    dynamic_data, n_latent = 0, kernel = "dynamic",
    residual_family = residual_family
  )
  model <- MCMedit(model, "B", c(1, 1), "Earnings_AR")
  model <- MCMedit(model, "B", c(1, 2), "Hours_lag_to_Earnings")
  model <- MCMedit(model, "B", c(2, 1), "Earnings_lag_to_Hours")
  MCMedit(model, "B", c(2, 2), "Hours_AR")
}
fit_dynamic_model <- function(
    model, seed, n_starts, rprop_iters, lbfgs_iters = 80L) {
  MCMfit(
    model, dynamic_data, compute_se = TRUE,
    optimizers = c("rprop", "lbfgs"),
    optim_iters = c(rprop_iters, lbfgs_iters),
    learning_rate = c(0.01, 0.005),
    moment_weighting = "diagonal", se_correction = "robust",
    n_starts = n_starts, seed = seed, verbose = FALSE
  )
}

message("Fitting dynamic MCMSEM with a Gaussian residual (20 starts)")
gaussian_model <- make_dynamic_model("gaussian")
gaussian_fit <- fit_dynamic_model(gaussian_model, 20260731L, 20L, 1400L)

constrain_signed_gamma_innovations <- function(model) {
  model <- MCMparameter(
    model, "shape_Earnings", "free", start = 0.114,
    transform = "positive"
  )
  model <- MCMparameter(
    model, "shape_Hours", "free", start = 0.295,
    transform = "positive"
  )
  model <- MCMparameter(model, "sign_Earnings", "fixed", value = -1)
  model <- MCMparameter(model, "sign_Hours", "fixed", value = 1)
  model <- MCMparameter(
    model, "tau_Earnings", "derived",
    expression = ~ sign_Earnings * 2 / sqrt(shape_Earnings)
  )
  model <- MCMparameter(
    model, "kappa_Earnings", "derived",
    expression = ~ 6 / shape_Earnings
  )
  model <- MCMparameter(
    model, "tau_Hours", "derived",
    expression = ~ sign_Hours * 2 / sqrt(shape_Hours)
  )
  MCMparameter(
    model, "kappa_Hours", "derived", expression = ~ 6 / shape_Hours
  )
}

gamma_grid <- data.frame(
  loading_earnings = c(1.32, 1.32, 1.32, -1.32, -1.32, -1.32, 0.25, -0.25),
  loading_hours = c(0.31, 0.31, 0.31, -0.31, -0.31, -0.31, 0.10, -0.10),
  shape = c(0.25, 2, 20, 0.25, 2, 20, 2, 2)
)
gamma_grid_fits <- lapply(seq_len(nrow(gamma_grid)), function(index) {
  message(
    "Fitting common-gamma grid ", index, "/", nrow(gamma_grid),
    " (3 starts)"
  )
  model <- constrain_signed_gamma_innovations(
    make_dynamic_model("common_gamma")
  )
  bounded_paths <- list(
    Earnings_AR = c(0.68, 0, 0.98),
    Hours_lag_to_Earnings = c(0.08, -0.98, 0.98),
    Earnings_lag_to_Hours = c(0.42, -0.98, 0.98),
    Hours_AR = c(0.65, 0, 0.98)
  )
  for (name in names(bounded_paths)) {
    specification <- bounded_paths[[name]]
    model <- MCMparameter(
      model, name, "free", start = specification[1L],
      transform = "bounded", lower = specification[2L],
      upper = specification[3L], overwrite = TRUE
    )
  }
  starts <- c(
    loading_Gamma_Earnings = gamma_grid$loading_earnings[index],
    loading_Gamma_Hours = gamma_grid$loading_hours[index],
    shape_Gamma = gamma_grid$shape[index]
  )
  for (name in names(starts)) {
    model <- MCMedit(model, "start", name, starts[[name]])
  }
  MCMfit(
    model, dynamic_data, compute_se = FALSE,
    optimizers = c("rprop", "lbfgs"), optim_iters = c(800, 60),
    learning_rate = c(0.01, 0.005),
    moment_weighting = "diagonal", se_correction = "robust",
    n_starts = 3L, seed = 20260830L + index, verbose = FALSE
  )
})
gamma_losses <- vapply(gamma_grid_fits, function(fit) fit$loss, numeric(1L))
gamma_best_grid <- which.min(gamma_losses)
gamma_fit <- MCMfit(
  gamma_grid_fits[[gamma_best_grid]], dynamic_data, compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"), optim_iters = c(500, 100),
  learning_rate = c(0.005, 0.002),
  moment_weighting = "diagonal", se_correction = "robust",
  n_starts = 1L, seed = 20260850L, verbose = FALSE
)

unrestricted_gamma_fit <- NULL
if (run_unrestricted_gamma) {
  unrestricted_grid <- expand.grid(
    orientation = c(-1, 1), shape = c(0.25, 2, 20, 150)
  )
  unrestricted_fits <- lapply(
    seq_len(nrow(unrestricted_grid)), function(index) {
      message(
        "Fitting unrestricted common-gamma grid ", index, "/",
        nrow(unrestricted_grid), " (3 starts)"
      )
      model <- make_dynamic_model("common_gamma")
      starts <- c(
        Earnings_AR = 0.676, Hours_lag_to_Earnings = 0.075,
        Earnings_lag_to_Hours = 0.422, Hours_AR = 0.651,
        loading_Gamma_Earnings =
          unrestricted_grid$orientation[index] * 1.35,
        loading_Gamma_Hours = unrestricted_grid$orientation[index] * 0.36,
        shape_Gamma = unrestricted_grid$shape[index]
      )
      for (name in names(starts)) {
        model <- MCMedit(model, "start", name, starts[[name]])
      }
      MCMfit(
        model, dynamic_data, compute_se = FALSE,
        optimizers = c("rprop", "lbfgs"), optim_iters = c(800, 60),
        learning_rate = c(0.01, 0.005),
        moment_weighting = "diagonal", se_correction = "robust",
        n_starts = 3L, seed = 20260902L + index, verbose = FALSE
      )
    }
  )
  unrestricted_losses <- vapply(
    unrestricted_fits, function(fit) fit$loss, numeric(1L)
  )
  unrestricted_best <- which.min(unrestricted_losses)
  unrestricted_gamma_fit <- MCMfit(
    unrestricted_fits[[unrestricted_best]], dynamic_data, compute_se = TRUE,
    optimizers = c("rprop", "lbfgs"), optim_iters = c(300, 80),
    learning_rate = c(0.005, 0.002),
    moment_weighting = "diagonal", se_correction = "robust",
    n_starts = 1L, seed = 20260920L, verbose = FALSE
  )
}

gaussian_paths <- extract_dynamic_paths(gaussian_fit)
gamma_paths <- extract_dynamic_paths(gamma_fit)
unrestricted_gamma_paths <- if (is.null(unrestricted_gamma_fit)) {
  NULL
} else {
  extract_dynamic_paths(unrestricted_gamma_fit)
}

fit_metrics <- data.frame(
  panel_rows = nrow(wide), selected_mcmsem_year = 2012L + mcm_wave,
  mcmsem_complete_n = nrow(mcm_values),
  clpm_cfi_scaled = lavaan::fitMeasures(clpm_fit, "cfi.scaled"),
  clpm_tli_scaled = lavaan::fitMeasures(clpm_fit, "tli.scaled"),
  clpm_rmsea_scaled = lavaan::fitMeasures(clpm_fit, "rmsea.scaled"),
  clpm_srmr = lavaan::fitMeasures(clpm_fit, "srmr"),
  riclpm_cfi_scaled = lavaan::fitMeasures(riclpm_fit, "cfi.scaled"),
  riclpm_tli_scaled = lavaan::fitMeasures(riclpm_fit, "tli.scaled"),
  riclpm_rmsea_scaled = lavaan::fitMeasures(riclpm_fit, "rmsea.scaled"),
  riclpm_srmr = lavaan::fitMeasures(riclpm_fit, "srmr"),
  gaussian_loss = gaussian_fit$loss,
  gaussian_information_condition = gaussian_fit$info$information_condition,
  gamma_loss = gamma_fit$loss,
  gamma_information_condition = gamma_fit$info$information_condition
)
if (!is.null(unrestricted_gamma_fit)) {
  fit_metrics$unrestricted_gamma_loss <- unrestricted_gamma_fit$loss
  fit_metrics$unrestricted_gamma_information_condition <-
    unrestricted_gamma_fit$info$information_condition
}

cat("Simulation truth\n")
print(simulation_truth)
cat("\nSimulation fit metrics\n")
print(simulation_metrics, row.names = FALSE)
cat("\nSimulation CLPM transition paths\n")
print(simulation_clpm_paths, row.names = FALSE)
cat("\nSimulation dynamic MCMSEM transition paths\n")
print(simulation_dynamic_paths, row.names = FALSE)
cat("\nGaussian-confounder simulation fit metrics\n")
print(simulation_gaussian_metrics, row.names = FALSE)
cat("\nGaussian-confounder simulation CLPM transition paths\n")
print(simulation_confounded_clpm_paths, row.names = FALSE)
cat("\nGaussian-confounder simulation RI-CLPM transition paths\n")
print(simulation_riclpm_paths, row.names = FALSE)
cat("\nGaussian-confounder simulation MCMSEM transition paths\n")
print(simulation_gaussian_paths, row.names = FALSE)
cat("\nGaussian-confounder simulation residual covariance\n")
print(simulation_gaussian_fit$Psi_G)
cat("\nReal-data fit metrics\n")
print(fit_metrics, row.names = FALSE)
cat("\nReal-data wave diagnostics\n")
print(wave_diagnostics, row.names = FALSE)
cat("\nCLPM transition paths\n")
print(clpm_paths, row.names = FALSE)
cat("\nRI-CLPM within-person transition paths\n")
print(riclpm_paths, row.names = FALSE)
cat("\nDynamic MCMSEM paths with Gaussian residual\n")
print(gaussian_paths, row.names = FALSE)
cat("\nEstimated Gaussian residual covariance\n")
print(gaussian_fit$Psi_G)
cat("\nDynamic MCMSEM paths with common-gamma residual\n")
print(gamma_paths, row.names = FALSE)
cat("\nEstimated common-gamma residual\n")
print(gamma_fit$dynamic$common_gamma)
if (!is.null(unrestricted_gamma_fit)) {
  cat("\nUnrestricted-innovation common-gamma sensitivity paths\n")
  print(unrestricted_gamma_paths, row.names = FALSE)
  cat("\nUnrestricted-innovation common-gamma residual\n")
  print(unrestricted_gamma_fit$dynamic$common_gamma)
}
output_dir <- file.path("validation-output", "longitudinal-example")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  as.data.frame(as.table(simulation_truth)),
  file.path(output_dir, "simulation_truth.csv"), row.names = FALSE
)
utils::write.csv(
  simulation_metrics,
  file.path(output_dir, "simulation_fit_metrics.csv"), row.names = FALSE
)
utils::write.csv(
  simulation_clpm_paths,
  file.path(output_dir, "simulation_clpm_paths.csv"), row.names = FALSE
)
utils::write.csv(
  simulation_dynamic_paths,
  file.path(output_dir, "simulation_dynamic_paths.csv"), row.names = FALSE
)
utils::write.csv(
  simulation_gaussian_metrics,
  file.path(output_dir, "simulation_gaussian_fit_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  simulation_confounded_clpm_paths,
  file.path(output_dir, "simulation_confounded_clpm_paths.csv"),
  row.names = FALSE
)
utils::write.csv(
  simulation_riclpm_paths,
  file.path(output_dir, "simulation_riclpm_paths.csv"), row.names = FALSE
)
utils::write.csv(
  simulation_gaussian_paths,
  file.path(output_dir, "simulation_gaussian_dynamic_paths.csv"),
  row.names = FALSE
)
utils::write.csv(
  simulation_gaussian_truth,
  file.path(output_dir, "simulation_gaussian_residual_truth.csv"),
  row.names = TRUE
)
utils::write.csv(
  simulation_gaussian_fit$Psi_G,
  file.path(output_dir, "simulation_gaussian_residual_estimate.csv"),
  row.names = TRUE
)
utils::write.csv(fit_metrics, file.path(output_dir, "fit_metrics.csv"), row.names = FALSE)
utils::write.csv(
  wave_diagnostics,
  file.path(output_dir, "wave_diagnostics.csv"), row.names = FALSE
)
utils::write.csv(clpm_paths, file.path(output_dir, "clpm_paths.csv"), row.names = FALSE)
utils::write.csv(
  riclpm_paths, file.path(output_dir, "riclpm_paths.csv"), row.names = FALSE
)
utils::write.csv(
  gaussian_paths,
  file.path(output_dir, "gaussian_dynamic_paths.csv"), row.names = FALSE
)
utils::write.csv(
  gaussian_fit$Psi_G,
  file.path(output_dir, "gaussian_residual_covariance.csv"), row.names = TRUE
)
utils::write.csv(
  gaussian_fit$start_diagnostics,
  file.path(output_dir, "gaussian_dynamic_start_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  gamma_paths,
  file.path(output_dir, "common_gamma_dynamic_paths.csv"), row.names = FALSE
)
utils::write.csv(
  gamma_fit$start_diagnostics,
  file.path(output_dir, "common_gamma_dynamic_start_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  transform(gamma_grid, loss = gamma_losses),
  file.path(output_dir, "common_gamma_start_grid.csv"), row.names = FALSE
)
if (!is.null(unrestricted_gamma_fit)) {
  utils::write.csv(
    unrestricted_gamma_paths,
    file.path(output_dir, "unrestricted_common_gamma_paths.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    transform(unrestricted_grid, loss = unrestricted_losses),
    file.path(output_dir, "unrestricted_common_gamma_start_grid.csv"),
    row.names = FALSE
  )
}
