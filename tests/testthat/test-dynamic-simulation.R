test_that("fast cross-sectional simulation recovers all dynamic parameter groups", {
  B_true <- matrix(c(0.40, 0.12, -0.08, 0.30), 2, byrow = TRUE)
  Psi_true <- matrix(c(0.25, 0.08, 0.08, 0.20), 2)
  data <- simulate_dynamic_cross_section(
    5000, B_true, Psi_G = Psi_true, seed = 88
  )
  summary_data <- MCMdatasummary(
    data, scale_data = FALSE, prep_asymptotic_se = FALSE
  )
  model <- MCMmodel(summary_data, kernel = "dynamic")
  L_true <- MCMSEM:::.covariance_to_cholesky_parameters(Psi_true)
  truth <- c(
    phi_X = 0.40, X_lag_to_Y = -0.08,
    Y_lag_to_X = 0.12, phi_Y = 0.30,
    tau_X = 1, tau_Y = -sqrt(8 / 5),
    kappa_X = 1.5, kappa_Y = 2.4,
    log_sd_G_X = L_true[1, 1], chol_G_Y_X = L_true[2, 1],
    log_sd_G_Y = L_true[2, 2]
  )
  true_loss <- MCMdynamicloss(model, summary_data, truth)
  fit <- MCMfit(
    model, summary_data, compute_se = FALSE,
    optimizers = c("rprop", "lbfgs"), optim_iters = c(150, 8),
    learning_rate = c(0.01, 0.2), n_starts = 2, seed = 88
  )

  expect_lt(max(abs(fit$B - B_true)), 0.20)
  expect_lt(max(abs(fit$innovation_third - c(X = 1, Y = -sqrt(8 / 5)))), 0.30)
  expect_lt(max(abs(fit$innovation_fourth - c(X = 1.5, Y = 2.4))), 1.0)
  expect_lt(max(abs(fit$Psi_G - Psi_true)), 0.25)
  expect_true(is.finite(true_loss))
  expect_lt(fit$loss, true_loss)
  expect_equal(nrow(fit$start_diagnostics), 2)
  expect_true(all(fit$start_diagnostics$stationary))
})
