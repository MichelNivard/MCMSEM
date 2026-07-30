test_that("dynamic residual-family API is explicit and backward compatible", {
  data <- make_moment_summary()
  legacy <- MCMmodel(data, kernel = "dynamic")
  none <- MCMmodel(data, kernel = "dynamic", gaussian_residual = FALSE)
  gamma <- MCMmodel(data, kernel = "dynamic", residual_family = "common_gamma")
  gamma_alias <- MCMmodel(data, kernel = "dynamic", residual_family = "gamma")

  expect_identical(MCMSEM:::.dynamic_residual_family(legacy), "gaussian")
  expect_identical(MCMSEM:::.dynamic_residual_family(none), "none")
  expect_identical(MCMSEM:::.dynamic_residual_family(gamma), "common_gamma")
  expect_identical(MCMSEM:::.dynamic_residual_family(gamma_alias), "common_gamma")
  expect_true(legacy$meta_data$gaussian_residual)
  expect_false(gamma$meta_data$gaussian_residual)
  expect_named(gamma$num_matrices, c(
    "B", "Tau", "Kappa", "L_G", "Lambda_Gamma", "Shape_Gamma", "D2"
  ))
  parameters <- MCMparameters(gamma)
  expect_equal(
    parameters$type[parameters$name == "shape_Gamma"], "free"
  )
  expect_equal(
    parameters$transform[parameters$name == "shape_Gamma"], "positive"
  )
  expect_equal(c(
    MCMdegreesoffreedom(legacy)$n_parameters,
    MCMdegreesoffreedom(gamma)$n_parameters,
    MCMdegreesoffreedom(gamma)$df
  ), c(11, 11, 1))

  expect_error(
    MCMmodel(
      data, kernel = "dynamic", gaussian_residual = TRUE,
      residual_family = "common_gamma"
    ),
    "conflicts"
  )
  expect_error(
    MCMmodel(data, residual_family = "common_gamma"),
    "only for.*dynamic"
  )
  expect_error(
    MCMmodel(data, kernel = "dynamic", residual_family = "lognormal"),
    "Invalid.*residual_family"
  )
})

test_that("common gamma contributes the analytic rank-one cumulants", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  values <- c(
    phi_A = 0, A_lag_to_B = 0, B_lag_to_A = 0, phi_B = 0,
    tau_A = 0, tau_B = 0, kappa_A = 0, kappa_B = 0,
    loading_Gamma_A = 0.6, loading_Gamma_B = -0.3,
    shape_Gamma = 4
  )
  model <- set_dynamic_parameters(model, values)
  implied <- MCMimpliedmoments(model)
  lambda <- c(0.6, -0.3)
  expected_M2 <- diag(2) + tcrossprod(lambda)
  expected_M3 <- MCMSEM:::.dynamic_array_to_mcm(
    MCMSEM:::.dynamic_vector_outer_power(lambda, 3L)
  )
  expected_K4 <- MCMSEM:::.dynamic_array_to_mcm(
    1.5 * MCMSEM:::.dynamic_vector_outer_power(lambda, 4L)
  )

  expect_equal(implied$M2, expected_M2, tolerance = 1e-12)
  expect_equal(implied$M3, expected_M3, tolerance = 1e-12)
  expect_equal(implied$K4, expected_K4, tolerance = 1e-12)
  expect_equal(unname(implied$residual_M2), tcrossprod(lambda), tolerance = 1e-12)
  expect_equal(implied$residual_M3, expected_M3, tolerance = 1e-12)
  expect_equal(implied$residual_K4, expected_K4, tolerance = 1e-12)
  expect_equal(unname(implied$Psi_G), matrix(0, 2, 2))
  expect_equal(implied$residual_shape, 4)
  expect_equal(implied$residual_skewness, 1)
  expect_equal(implied$residual_excess_kurtosis, 1.5)
  expect_equal(implied$residual_loadings, c(A = 0.6, B = -0.3))
  expect_equal(
    implied$M4,
    MCMSEM:::.dynamic_array_to_mcm(
      MCMSEM:::.cumulant4_to_raw4(
        MCMSEM:::.dynamic_mcm_to_array(expected_K4, 4L), expected_M2
      )
    ),
    tolerance = 1e-12
  )
})

test_that("common-gamma base and Torch moments and gradients agree", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  model <- set_dynamic_parameters(model, c(
    phi_A = 0.3, A_lag_to_B = -0.1, B_lag_to_A = 0.15,
    phi_B = 0.25, tau_A = 0.7, tau_B = -0.5,
    kappa_A = 2, kappa_B = 3,
    loading_Gamma_A = 0.6, loading_Gamma_B = -0.3,
    shape_Gamma = 4
  ))
  components <- MCMSEM:::.get_dynamic_torch_components(
    model, torch_device("cpu"), torch_float64()
  )
  torch_implied <- MCMSEM:::.get_dynamic_predicted_matrices(components)
  base_implied <- MCMimpliedmoments(model)
  for (name in c(
    "M2", "M3", "M4", "K4", "residual_M2", "residual_M3", "residual_K4"
  )) {
    expect_equal(
      as.matrix(torch_implied[[name]]), unname(base_implied[[name]]),
      tolerance = 1e-10
    )
  }
  objective <- torch_sum(torch_implied$M2) + torch_sum(torch_implied$M3) +
    torch_sum(torch_implied$M4)
  objective$backward()
  gradients <- unlist(MCMSEM:::.dynamic_gradient_snapshot(components))
  expect_true(all(is.finite(gradients)))
  expect_true(all(c(
    "graph.loading_Gamma_A", "graph.loading_Gamma_B", "graph.shape_Gamma"
  ) %in% names(gradients)))
  expect_true(all(abs(gradients[c(
    "graph.loading_Gamma_A", "graph.loading_Gamma_B", "graph.shape_Gamma"
  )]) > 1e-8))
})

test_that("common-gamma moments are permutation equivariant", {
  variables <- c("A", "B", "C")
  model <- MCMmodel(
    make_moment_summary(3, variables), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  B <- matrix(c(
    .35, .10, -.04, -.06, .30, .08, .03, -.05, .32
  ), 3, byrow = TRUE)
  model$num_matrices$B <- B
  model$num_matrices$Tau[, 1] <- c(-.8, .6, 1.1)
  model$num_matrices$Kappa[, 1] <- c(2, 3, 4)
  model$num_matrices$Lambda_Gamma[, 1] <- c(.5, -.2, .35)
  model$num_matrices$Shape_Gamma[1, 1] <- 2.5
  first <- MCMimpliedmoments(model)

  permutation <- c(3, 1, 2)
  other <- MCMmodel(
    make_moment_summary(3, variables[permutation]), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  other$num_matrices$B <- B[permutation, permutation]
  other$num_matrices$Tau[, 1] <- model$num_matrices$Tau[permutation, 1]
  other$num_matrices$Kappa[, 1] <- model$num_matrices$Kappa[permutation, 1]
  other$num_matrices$Lambda_Gamma[, 1] <-
    model$num_matrices$Lambda_Gamma[permutation, 1]
  other$num_matrices$Shape_Gamma[1, 1] <- 2.5
  second <- MCMimpliedmoments(other)

  expect_equal(second$M2, first$M2[permutation, permutation], tolerance = 1e-10)
  expect_equal(
    MCMSEM:::.dynamic_mcm_to_array(second$M3, 3L),
    MCMSEM:::.dynamic_mcm_to_array(first$M3, 3L)[
      permutation, permutation, permutation
    ], tolerance = 1e-10
  )
  expect_equal(
    MCMSEM:::.dynamic_mcm_to_array(second$M4, 4L),
    MCMSEM:::.dynamic_mcm_to_array(first$M4, 4L)[
      permutation, permutation, permutation, permutation
    ], tolerance = 1e-9
  )
})

test_that("common-gamma population fit reports components and robust SEs", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  model <- set_dynamic_parameters(model, c(
    phi_A = 0.45, A_lag_to_B = -0.12, B_lag_to_A = 0.18,
    phi_B = 0.35, tau_A = 1, tau_B = -0.8,
    kappa_A = 2, kappa_B = 4,
    loading_Gamma_A = 0.5, loading_Gamma_B = 0.25,
    shape_Gamma = 3
  ))
  population <- prepare_population_vcov(make_population_summary(model))
  fit <- MCMfit(
    model, population, compute_se = TRUE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1,
    n_starts = 1, moment_weighting = "diagonal", weight_ridge = 0
  )
  expect_lt(fit$loss, 1e-16)
  expect_identical(fit$info$residual_family, "common_gamma")
  expect_identical(fit$dynamic$residual_family, "common_gamma")
  expect_equal(
    unname(fit$dynamic$residual_M2), tcrossprod(c(.5, .25)), tolerance = 1e-8
  )
  expect_equal(fit$dynamic$common_gamma$loadings$estimate, c(.5, .25),
               tolerance = 1e-8)
  expect_equal(fit$dynamic$common_gamma$shape$estimate, 3, tolerance = 1e-8)
  expect_true(all(is.finite(
    unlist(fit$dynamic$common_gamma$shape[
      c("se", "skewness_se", "excess_kurtosis_se")
    ], use.names = FALSE)
  )))
  expect_true(all(is.finite(fit$dynamic$residual_covariance$se)))
  expect_equal(nrow(fit$dynamic$gaussian_covariance), 0L)
  expect_equal(fit$dynamic$asymptotic$jacobian_rank, 11)
  expect_true(fit$info$standard_errors_available)
  expect_output(print(summary(fit)), "Nominal df     : 1")
})

test_that("gamma innovations plus common-gamma confounder have three df", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic",
    residual_family = "common_gamma"
  )
  for (variable in c("A", "B")) {
    model <- MCMparameter(
      model, paste0("shape_", variable), "free", start = 4,
      transform = "positive"
    )
    model <- MCMparameter(
      model, paste0("sign_", variable), "fixed",
      value = if (variable == "A") -1 else 1
    )
  }
  model <- MCMparameter(
    model, "tau_A", "derived", expression = ~ sign_A * 2 / sqrt(shape_A)
  )
  model <- MCMparameter(
    model, "kappa_A", "derived", expression = ~ 6 / shape_A
  )
  model <- MCMparameter(
    model, "tau_B", "derived", expression = ~ sign_B * 2 / sqrt(shape_B)
  )
  model <- MCMparameter(
    model, "kappa_B", "derived", expression = ~ 6 / shape_B
  )
  degrees <- MCMdegreesoffreedom(model)
  expect_equal(c(degrees$n_parameters, degrees$df), c(9, 3))
  expect_equal(MCMdiagnostics(model, jacobian = TRUE)$jacobian_rank, 9)
})
