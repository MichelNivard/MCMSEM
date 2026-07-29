test_that("population dynamic fit uses the common result infrastructure", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  truth <- c(
    phi_A = 0.45, A_lag_to_B = -0.12, B_lag_to_A = 0.18,
    phi_B = 0.35, tau_A = 1, tau_B = -0.8,
    kappa_A = 2, kappa_B = 4,
    log_sd_G_A = log(0.5), chol_G_B_A = 0.2,
    log_sd_G_B = log(0.6)
  )
  model <- set_dynamic_parameters(model, truth)
  population <- make_population_summary(model)
  fit <- MCMfit(
    model, population, compute_se = FALSE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1,
    n_starts = 1, seed = 42
  )
  expect_s4_class(fit, "mcmresultclass")
  expect_identical(fit$kernel, "dynamic")
  expect_identical(fit$model$meta_data$kernel, "dynamic")
  expect_lt(fit$loss, 1e-20)
  expect_equal(fit$B, model$num_matrices$B, tolerance = 1e-8)
  expect_equal(unname(fit$Psi_G), MCMSEM:::.dynamic_gaussian_covariance(model$num_matrices$L_G),
               tolerance = 1e-8)
  expect_equal(fit$innovation_variances, c(A = 1, B = 1))
  expect_equal(c(fit$n_moments, fit$degrees_of_freedom), c(12, 1))
  expect_true(fit$stationary)
  expect_named(fit$predicted,
               c("M2", "M3", "M4", "K4", "within_M2", "Psi_G", "L_G"))
  expect_named(fit$residuals, c("M2", "M3", "M4"))
  expect_output(print(fit), "Kernel: Stationary Dynamic MCMSEM")
  expect_output(print(summary(fit)), "Nominal df     : 1")
})

test_that("dynamic robust WLS standard errors use the sandwich equation", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  truth <- c(
    phi_A = 0.45, A_lag_to_B = -0.12, B_lag_to_A = 0.18,
    phi_B = 0.35, tau_A = 1, tau_B = -0.8,
    kappa_A = 2, kappa_B = 4,
    log_sd_G_A = log(0.5), chol_G_B_A = 0.2,
    log_sd_G_B = log(0.6)
  )
  model <- set_dynamic_parameters(model, truth)
  population <- make_population_summary(model)
  A <- diag(12)
  A[row(A) != col(A)] <- 0.03
  Omega <- (A %*% t(A)) / population$meta_data$N
  population <- prepare_population_vcov(population, Omega)
  fit <- MCMfit(
    model, population, compute_se = TRUE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1,
    n_starts = 1, moment_weighting = "diagonal", weight_ridge = 0
  )
  ase <- fit$dynamic$asymptotic
  Delta <- ase$jacobian
  W <- fit$dynamic$moment_weight
  bread <- solve(crossprod(Delta, W %*% Delta))
  expected <- bread %*% crossprod(Delta, W %*% Omega %*% W %*% Delta) %*% bread
  expect_equal(unname(ase$vcov), unname(expected), tolerance = 1e-6)
  expect_identical(fit$info$se_correction, "robust")
  expect_true(fit$info$standard_errors_available)
  expect_true(all(is.finite(unlist(fit$df["se", ]))))
  expect_true(all(is.finite(fit$dynamic$gaussian_covariance$se)))
  expect_equal(ase$jacobian_rank, length(model$param_values))
})

test_that("full WLS uses the efficient information covariance", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- set_dynamic_parameters(model, c(
    phi_A = 0.45, A_lag_to_B = -0.12, B_lag_to_A = 0.18,
    phi_B = 0.35, tau_A = 1, tau_B = -0.8,
    kappa_A = 2, kappa_B = 4,
    log_sd_G_A = log(0.5), chol_G_B_A = 0.2,
    log_sd_G_B = log(0.6)
  ))
  population <- prepare_population_vcov(make_population_summary(model))
  specification <- MCMSEM:::.dynamic_weight_specification(
    population, "full", weight_ridge = 0, require_vcov = TRUE
  )
  ase <- MCMSEM:::.dynamic_asymptotic_se(
    model, population, specification, se_correction = "auto",
    weight_ridge = 0
  )
  expected <- solve(crossprod(ase$jacobian, solve(diag(12) / 10000) %*%
                                ase$jacobian))
  expect_identical(ase$correction, "model_based")
  expect_equal(unname(ase$vcov), unname(expected), tolerance = 1e-6)
})

test_that("dynamic SE rejects old or unprepared moment summaries", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  population <- make_population_summary(model)
  expect_error(
    MCMfit(
      model, population, compute_se = TRUE,
      optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1
    ),
    "prep_asymptotic_se"
  )
})

test_that("dynamic model without Gaussian residual reports no Gaussian parameters", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic", gaussian_residual = FALSE
  )
  population <- prepare_population_vcov(make_population_summary(model))
  fit <- MCMfit(
    model, population, compute_se = TRUE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1
  )
  expect_equal(nrow(fit$dynamic$gaussian_covariance), 0L)
  expect_identical(dim(fit$dynamic$asymptotic$gaussian_vcov), c(0L, 0L))
  expect_length(fit$dynamic$asymptotic$gaussian_se, 0L)
})

test_that("multi-start generation is seeded and stationary", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  first <- MCMSEM:::.dynamic_random_starts(model, 4, seed = 99)
  second <- MCMSEM:::.dynamic_random_starts(model, 4, seed = 99)
  expect_equal(first, second)
  for (start in first) {
    candidate <- MCMSEM:::.dynamic_model_with_parameters(model, start)
    expect_lt(MCMSEM:::.dynamic_spectral_radius(candidate$num_matrices$B), 0.995)
    expect_gte(candidate$num_matrices$B[1, 1], 0)
    expect_gte(candidate$num_matrices$B[2, 2], 0)
  }
})
