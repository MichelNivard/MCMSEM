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
simulation_dynamic_paths <- simulation_dynamic_fit$transition_parameters[
  , c("label", "lagged", "current", "estimate", "se")
]
simulation_dynamic_paths$pvalue <- 2 * stats::pnorm(
  abs(simulation_dynamic_paths$estimate / simulation_dynamic_paths$se),
  lower.tail = FALSE
)
simulation_dynamic_paths$ci.lower <- simulation_dynamic_paths$estimate -
  1.96 * simulation_dynamic_paths$se
simulation_dynamic_paths$ci.upper <- simulation_dynamic_paths$estimate +
  1.96 * simulation_dynamic_paths$se
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

# Then run the same comparison on real longitudinal data.
data_url <- paste0(
  "https://vincentarelbundock.github.io/Rdatasets/csv/",
  "sampleSelection/nlswork.csv"
)
data_file <- tempfile(fileext = ".csv")
on.exit(unlink(data_file), add = TRUE)
utils::download.file(data_url, data_file, mode = "wb", quiet = TRUE)
expected_md5 <- "f546ffe0bee86acb5d79b8775d341709"
observed_md5 <- unname(tools::md5sum(data_file))
if (!identical(observed_md5, expected_md5)) {
  stop("The downloaded NLS data do not match the validated file version.")
}

nls_long <- utils::read.csv(data_file)
nls_long <- nls_long[nls_long$year %in% c(71, 73, 75, 77), ]

# MCMSEM fixes both innovation variances to one. These transparent linear
# rescalings retain all distributional information while putting both observed
# variances above one. They are used in the CLPM and MCMSEM fits alike.
nls_long$Wage <- 4 * nls_long$ln_wage
nls_long$Hours <- nls_long$hours / 5

wide <- stats::reshape(
  nls_long[c("idcode", "year", "Wage", "Hours")],
  idvar = "idcode", timevar = "year", direction = "wide"
)
names(wide) <- sub("Wage\\.", "Wage", names(wide))
names(wide) <- sub("Hours\\.", "Hours", names(wide))

clpm_syntax <- "
  Wage73 ~ wage_ar*Wage71 + hours_to_wage*Hours71
  Wage75 ~ wage_ar*Wage73 + hours_to_wage*Hours73
  Wage77 ~ wage_ar*Wage75 + hours_to_wage*Hours75

  Hours73 ~ hours_ar*Hours71 + wage_to_hours*Wage71
  Hours75 ~ hours_ar*Hours73 + wage_to_hours*Wage73
  Hours77 ~ hours_ar*Hours75 + wage_to_hours*Wage75

  Wage71 ~~ Hours71
  Wage73 ~~ Hours73
  Wage75 ~~ Hours75
  Wage77 ~~ Hours77
"

clpm_fit <- lavaan::sem(
  clpm_syntax, data = wide,
  estimator = "MLR", missing = "fiml", meanstructure = TRUE
)
clpm_estimates <- lavaan::parameterEstimates(
  clpm_fit, standardized = TRUE, ci = TRUE
)
clpm_labels <- c("wage_ar", "hours_to_wage", "hours_ar", "wage_to_hours")
clpm_paths <- clpm_estimates[
  !duplicated(clpm_estimates$label) & clpm_estimates$label %in% clpm_labels,
  c("label", "est", "se", "pvalue", "std.all", "ci.lower", "ci.upper")
]
clpm_paths <- clpm_paths[match(clpm_labels, clpm_paths$label), ]

riclpm_syntax <- "
  RI_Wage =~ 1*Wage71 + 1*Wage73 + 1*Wage75 + 1*Wage77
  RI_Hours =~ 1*Hours71 + 1*Hours73 + 1*Hours75 + 1*Hours77

  wWage71 =~ 1*Wage71
  wWage73 =~ 1*Wage73
  wWage75 =~ 1*Wage75
  wWage77 =~ 1*Wage77
  wHours71 =~ 1*Hours71
  wHours73 =~ 1*Hours73
  wHours75 =~ 1*Hours75
  wHours77 =~ 1*Hours77

  Wage71 ~~ 0*Wage71
  Wage73 ~~ 0*Wage73
  Wage75 ~~ 0*Wage75
  Wage77 ~~ 0*Wage77
  Hours71 ~~ 0*Hours71
  Hours73 ~~ 0*Hours73
  Hours75 ~~ 0*Hours75
  Hours77 ~~ 0*Hours77

  wWage73 ~ wage_ar*wWage71 + hours_to_wage*wHours71
  wWage75 ~ wage_ar*wWage73 + hours_to_wage*wHours73
  wWage77 ~ wage_ar*wWage75 + hours_to_wage*wHours75
  wHours73 ~ hours_ar*wHours71 + wage_to_hours*wWage71
  wHours75 ~ hours_ar*wHours73 + wage_to_hours*wWage73
  wHours77 ~ hours_ar*wHours75 + wage_to_hours*wWage75

  RI_Wage ~~ RI_Hours
  wWage71 ~~ wHours71
  wWage73 ~~ wHours73
  wWage75 ~~ wHours75
  wWage77 ~~ wHours77
  RI_Wage ~~ 0*wWage71 + 0*wHours71
  RI_Hours ~~ 0*wWage71 + 0*wHours71
"
riclpm_fit <- lavaan::sem(
  riclpm_syntax, data = wide,
  estimator = "MLR", missing = "fiml", meanstructure = TRUE,
  fixed.x = FALSE
)
if (!isTRUE(lavaan::lavInspect(riclpm_fit, "post.check"))) {
  stop("The RI-CLPM solution failed lavaan's post-estimation check.")
}
riclpm_estimates <- lavaan::parameterEstimates(
  riclpm_fit, standardized = TRUE, ci = TRUE
)
riclpm_paths <- riclpm_estimates[
  !duplicated(riclpm_estimates$label) &
    riclpm_estimates$label %in% clpm_labels,
  c("label", "est", "se", "pvalue", "std.all", "ci.lower", "ci.upper")
]
riclpm_paths <- riclpm_paths[match(clpm_labels, riclpm_paths$label), ]

wave_years <- c(71, 73, 75, 77)
complete_n <- vapply(wave_years, function(yy) {
  sum(stats::complete.cases(wide[paste0(c("Wage", "Hours"), yy)]))
}, integer(1))
mcm_wave <- wave_years[which.max(complete_n)]
final_wave <- stats::na.omit(
  wide[paste0(c("Wage", "Hours"), mcm_wave)]
)
names(final_wave) <- c("Wage", "Hours")
wave_diagnostics <- do.call(rbind, lapply(wave_years, function(yy) {
  values <- stats::na.omit(wide[paste0(c("Wage", "Hours"), yy)])
  names(values) <- c("Wage", "Hours")
  data.frame(
    year = 1900 + yy,
    complete_n = nrow(values),
    wage_mean = mean(values$Wage),
    hours_mean = mean(values$Hours),
    wage_variance = stats::var(values$Wage),
    hours_variance = stats::var(values$Hours),
    covariance = stats::cov(values$Wage, values$Hours)
  )
}))
dynamic_data <- MCMdatasummary(
  final_wave,
  scale_data = FALSE,
  prep_asymptotic_se = TRUE,
  use_skewness = TRUE,
  use_kurtosis = TRUE
)
dynamic_model <- MCMmodel(
  dynamic_data,
  n_latent = 0,
  kernel = "dynamic",
  gaussian_residual = FALSE
)
dynamic_model <- MCMedit(dynamic_model, "B", c(1, 1), "Wage_AR")
dynamic_model <- MCMedit(dynamic_model, "B", c(1, 2), "Hours_lag_to_Wage")
dynamic_model <- MCMedit(dynamic_model, "B", c(2, 1), "Wage_lag_to_Hours")
dynamic_model <- MCMedit(dynamic_model, "B", c(2, 2), "Hours_AR")

dynamic_fit <- MCMfit(
  dynamic_model,
  dynamic_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(750, 40),
  learning_rate = c(0.01, 0.2),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 20,
  seed = 20260729,
  verbose = FALSE
)

gaussian_model <- MCMmodel(
  dynamic_data,
  n_latent = 0,
  kernel = "dynamic",
  gaussian_residual = TRUE
)
gaussian_model <- MCMedit(gaussian_model, "B", c(1, 1), "Wage_AR")
gaussian_model <- MCMedit(
  gaussian_model, "B", c(1, 2), "Hours_lag_to_Wage"
)
gaussian_model <- MCMedit(
  gaussian_model, "B", c(2, 1), "Wage_lag_to_Hours"
)
gaussian_model <- MCMedit(gaussian_model, "B", c(2, 2), "Hours_AR")
gaussian_fit <- MCMfit(
  gaussian_model,
  dynamic_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(1000, 50),
  learning_rate = c(0.01, 0.2),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 30,
  seed = 20260731,
  verbose = FALSE
)

dynamic_paths <- dynamic_fit$transition_parameters[
  , c("label", "lagged", "current", "estimate", "se")
]
dynamic_paths$pvalue <- 2 * stats::pnorm(
  abs(dynamic_paths$estimate / dynamic_paths$se), lower.tail = FALSE
)
dynamic_paths$ci.lower <- dynamic_paths$estimate - 1.96 * dynamic_paths$se
dynamic_paths$ci.upper <- dynamic_paths$estimate + 1.96 * dynamic_paths$se
gaussian_paths <- gaussian_fit$transition_parameters[
  , c("label", "lagged", "current", "estimate", "se")
]
gaussian_paths$pvalue <- 2 * stats::pnorm(
  abs(gaussian_paths$estimate / gaussian_paths$se), lower.tail = FALSE
)
gaussian_paths$ci.lower <- gaussian_paths$estimate -
  1.96 * gaussian_paths$se
gaussian_paths$ci.upper <- gaussian_paths$estimate +
  1.96 * gaussian_paths$se

fit_metrics <- data.frame(
  panel_participants = nrow(wide),
  selected_mcmsem_year = 1900 + mcm_wave,
  final_wave_complete = nrow(final_wave),
  clpm_cfi_robust = lavaan::fitMeasures(clpm_fit, "cfi.robust"),
  clpm_tli_robust = lavaan::fitMeasures(clpm_fit, "tli.robust"),
  clpm_rmsea_robust = lavaan::fitMeasures(clpm_fit, "rmsea.robust"),
  clpm_srmr = lavaan::fitMeasures(clpm_fit, "srmr"),
  riclpm_cfi_robust = lavaan::fitMeasures(riclpm_fit, "cfi.robust"),
  riclpm_tli_robust = lavaan::fitMeasures(riclpm_fit, "tli.robust"),
  riclpm_rmsea_robust = lavaan::fitMeasures(riclpm_fit, "rmsea.robust"),
  riclpm_srmr = lavaan::fitMeasures(riclpm_fit, "srmr"),
  dynamic_loss = dynamic_fit$loss,
  dynamic_spectral_radius = dynamic_fit$spectral_radius,
  dynamic_nominal_df = dynamic_fit$degrees_of_freedom,
  dynamic_jacobian_rank = dynamic_fit$info$jacobian_rank,
  dynamic_information_condition = dynamic_fit$info$information_condition,
  dynamic_admissible_starts = sum(dynamic_fit$start_diagnostics$convergence == 0),
  gaussian_dynamic_loss = gaussian_fit$loss,
  gaussian_dynamic_spectral_radius = gaussian_fit$spectral_radius,
  gaussian_dynamic_nominal_df = gaussian_fit$degrees_of_freedom,
  gaussian_dynamic_jacobian_rank = gaussian_fit$info$jacobian_rank,
  gaussian_dynamic_information_condition =
    gaussian_fit$info$information_condition,
  gaussian_dynamic_admissible_starts = sum(
    gaussian_fit$start_diagnostics$convergence == 0
  )
)

cat("Simulation truth\n")
print(simulation_truth)
cat("\nSimulation fit metrics\n")
print(simulation_metrics, row.names = FALSE)
cat("\nSimulation CLPM transition paths\n")
print(simulation_clpm_paths, row.names = FALSE)
cat("\nSimulation dynamic MCMSEM transition paths\n")
print(simulation_dynamic_paths, row.names = FALSE)
cat("\nReal-data fit metrics\n")
print(fit_metrics, row.names = FALSE)
cat("\nReal-data wave diagnostics\n")
print(wave_diagnostics, row.names = FALSE)
cat("\nCLPM transition paths\n")
print(clpm_paths, row.names = FALSE)
cat("\nRI-CLPM within-person transition paths\n")
print(riclpm_paths, row.names = FALSE)
cat("\nDynamic MCMSEM transition paths\n")
print(dynamic_paths, row.names = FALSE)
cat("\nDynamic MCMSEM paths with Gaussian residual\n")
print(gaussian_paths, row.names = FALSE)
cat("\nEstimated Gaussian residual covariance\n")
print(gaussian_fit$Psi_G)
cat("\nDynamic start diagnostics\n")
print(dynamic_fit$start_diagnostics, row.names = FALSE)

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
utils::write.csv(fit_metrics, file.path(output_dir, "fit_metrics.csv"), row.names = FALSE)
utils::write.csv(
  wave_diagnostics,
  file.path(output_dir, "wave_diagnostics.csv"), row.names = FALSE
)
utils::write.csv(clpm_paths, file.path(output_dir, "clpm_paths.csv"), row.names = FALSE)
utils::write.csv(
  riclpm_paths, file.path(output_dir, "riclpm_paths.csv"), row.names = FALSE
)
utils::write.csv(dynamic_paths, file.path(output_dir, "dynamic_paths.csv"), row.names = FALSE)
utils::write.csv(
  gaussian_paths,
  file.path(output_dir, "gaussian_dynamic_paths.csv"), row.names = FALSE
)
utils::write.csv(
  gaussian_fit$Psi_G,
  file.path(output_dir, "gaussian_residual_covariance.csv"), row.names = TRUE
)
utils::write.csv(
  dynamic_fit$start_diagnostics,
  file.path(output_dir, "dynamic_start_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  gaussian_fit$start_diagnostics,
  file.path(output_dir, "gaussian_dynamic_start_diagnostics.csv"),
  row.names = FALSE
)
