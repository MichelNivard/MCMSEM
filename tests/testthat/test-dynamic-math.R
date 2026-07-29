test_that("zero transition returns innovation cumulants", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic",
                    gaussian_residual = FALSE)
  values <- c(phi_A = 0, A_lag_to_B = 0, B_lag_to_A = 0, phi_B = 0,
              tau_A = 0.8, tau_B = -0.5, kappa_A = 2, kappa_B = 4)
  model <- set_dynamic_parameters(model, values)
  implied <- MCMSEM:::.dynamic_implied_moments_base(model)
  expect_equal(implied$M2, diag(2))
  expect_equal(diag(implied$M3[, c(1, 4)]), c(0.8, -0.5))
  expect_equal(c(implied$K4[1, 1], implied$K4[2, 8]), c(2, 4))
  expect_equal(unname(implied$Psi_G), matrix(0, 2, 2))
})

test_that("univariate cumulants have the analytic stationary form", {
  phi <- 0.6
  for (order in 2:4) {
    D <- array(2.5, dim = rep(1, order))
    actual <- MCMSEM:::.dynamic_cumulant_propagation(matrix(phi), D, order)
    expect_equal(actual, 2.5 / (1 - phi^order), tolerance = 1e-12)
  }
})

test_that("linear solve agrees with a long finite series", {
  B <- matrix(c(0.45, 0.12, -0.08, 0.35), 2, byrow = TRUE)
  for (order in 2:4) {
    D <- array(0, rep(2, order))
    D[matrix(rep(1, order), nrow = 1)] <- 1.2
    D[matrix(rep(2, order), nrow = 1)] <- -0.4
    K <- MCMSEM:::.dynamic_kron_power(B, order)
    term <- as.vector(D)
    finite <- term
    for (h in seq_len(250L)) {
      term <- K %*% term
      finite <- finite + term
    }
    solved <- MCMSEM:::.dynamic_cumulant_propagation(B, D, order)
    expect_equal(as.numeric(solved), as.numeric(finite), tolerance = 1e-11)
  }
})

test_that("Gaussian residuals affect covariance and raw fourth moments only", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  values <- c(phi_A = 0.4, A_lag_to_B = -0.1, B_lag_to_A = 0.15,
              phi_B = 0.3, tau_A = 1, tau_B = -0.7,
              kappa_A = 2, kappa_B = 3,
              log_sd_G_A = log(0.3), chol_G_B_A = 0.1,
              log_sd_G_B = log(0.4))
  model <- set_dynamic_parameters(model, values)
  first <- MCMSEM:::.dynamic_implied_moments_base(model)
  model <- MCMedit(model, "start", "log_sd_G_A", log(0.8))
  model <- MCMedit(model, "start", "chol_G_B_A", -0.2)
  second <- MCMSEM:::.dynamic_implied_moments_base(model)
  expect_false(isTRUE(all.equal(first$M2, second$M2)))
  expect_false(isTRUE(all.equal(first$M4, second$M4)))
  expect_equal(first$M3, second$M3, tolerance = 1e-12)
  expect_equal(first$K4, second$K4, tolerance = 1e-12)
})

test_that("raw fourth moments and cumulants round trip", {
  Sigma <- matrix(c(2, 0.4, 0.4, 1.3), 2)
  K4 <- array(seq_len(16) / 10, c(2, 2, 2, 2))
  raw <- MCMSEM:::.cumulant4_to_raw4(K4, Sigma)
  expect_equal(MCMSEM:::.raw4_to_cumulant4(raw, Sigma), K4,
               tolerance = 1e-12)
})

test_that("Cholesky parameters always produce symmetric PSD covariance", {
  set.seed(20)
  for (i in seq_len(20)) {
    L <- matrix(rnorm(16), 4)
    Psi <- MCMSEM:::.dynamic_gaussian_covariance(L)
    expect_equal(Psi, t(Psi), tolerance = 1e-12)
    expect_true(min(eigen(Psi, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
  }
})

test_that("stationarity is enforced and starts are screened", {
  expect_error(
    MCMSEM:::.dynamic_cumulant_propagation(diag(1.01, 2), diag(2), 2),
    "stationary"
  )
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMedit(model, "start", "phi_A", 0.98)
  model <- MCMedit(model, "start", "phi_B", 0.98)
  model <- MCMedit(model, "start", "A_lag_to_B", 0.98)
  model <- MCMedit(model, "start", "B_lag_to_A", 0.98)
  screened <- MCMSEM:::.dynamic_screen_start(model, model$param_values, 0.995)
  screened_model <- MCMSEM:::.dynamic_model_with_parameters(model, screened)
  expect_lt(MCMSEM:::.dynamic_spectral_radius(screened_model$num_matrices$B), 0.995)
})

test_that("dynamic moments are permutation equivariant", {
  model <- MCMmodel(make_moment_summary(3), kernel = "dynamic")
  B <- matrix(c(0.4, 0.1, -0.05, -0.08, 0.3, 0.12,
                0.04, -0.06, 0.35), 3, byrow = TRUE)
  tau <- c(0.7, -1.1, 0.4)
  kappa <- c(1.5, 3, -0.5)
  L <- matrix(c(log(.4), .1, -.2, 0, log(.5), .15, 0, 0, log(.6)), 3)
  model$num_matrices$B <- B
  model$num_matrices$Tau[, 1] <- tau
  model$num_matrices$Kappa[, 1] <- kappa
  model$num_matrices$L_G <- L
  first <- MCMSEM:::.dynamic_implied_moments_base(model)

  perm <- c(3, 1, 2)
  second_model <- MCMmodel(make_moment_summary(3, LETTERS[perm]), kernel = "dynamic")
  second_model$num_matrices$B <- B[perm, perm]
  second_model$num_matrices$Tau[, 1] <- tau[perm]
  second_model$num_matrices$Kappa[, 1] <- kappa[perm]
  Psi_perm <- first$Psi_G[perm, perm]
  second_model$num_matrices$L_G <-
    MCMSEM:::.covariance_to_cholesky_parameters(Psi_perm)
  second <- MCMSEM:::.dynamic_implied_moments_base(second_model)

  expect_equal(second$M2, first$M2[perm, perm], tolerance = 1e-10)
  expect_equal(
    MCMSEM:::.dynamic_mcm_to_array(second$M3, 3L),
    MCMSEM:::.dynamic_mcm_to_array(first$M3, 3L)[perm, perm, perm],
    tolerance = 1e-10
  )
  expect_equal(
    MCMSEM:::.dynamic_mcm_to_array(second$M4, 4L),
    MCMSEM:::.dynamic_mcm_to_array(first$M4, 4L)[perm, perm, perm, perm],
    tolerance = 1e-9
  )
})

test_that("population loss, moment count, and nominal df are correct", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  population <- make_population_summary(model)
  expect_lt(MCMdynamicloss(model, population), 1e-24)
  counts <- MCMmomentcount(2)
  expect_equal(unname(counts), c(3, 4, 5, 12))
  dof <- MCMdegreesoffreedom(model)
  expect_equal(c(dof$n_moments, dof$n_parameters, dof$df), c(12, 11, 1))
})

test_that("torch and base implementations agree and torch gradients are finite", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  components <- MCMSEM:::.get_dynamic_torch_components(
    model, torch_device("cpu"), torch_float64()
  )
  torch_implied <- MCMSEM:::.get_dynamic_predicted_matrices(components)
  base_implied <- MCMSEM:::.dynamic_implied_moments_base(model)
  expect_equal(as.matrix(torch_implied$M2), base_implied$M2, tolerance = 1e-11)
  expect_equal(as.matrix(torch_implied$M3), base_implied$M3, tolerance = 1e-11)
  expect_equal(as.matrix(torch_implied$M4), base_implied$M4, tolerance = 1e-10)
  objective <- torch_sum(torch_implied$M2) + torch_sum(torch_implied$M3) +
    torch_sum(torch_implied$M4)
  objective$backward()
  gradients <- unlist(MCMSEM:::.dynamic_gradient_snapshot(components))
  expect_true(all(is.finite(gradients)))
  expect_true(any(abs(gradients) > 1e-8))
})
