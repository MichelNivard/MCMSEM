test_that("kernel API is exact, canonical, and backward compatible", {
  data <- make_moment_summary()
  default_model <- MCMmodel(data)
  canonical_model <- MCMmodel(data, kernel = "contemporaneous")
  alias_model <- NULL
  expect_silent(alias_model <- MCMmodel(data, kernel = "static"))

  expect_identical(default_model$meta_data$kernel, "contemporaneous")
  expect_identical(canonical_model$meta_data$kernel, "contemporaneous")
  expect_identical(alias_model$meta_data$kernel, "contemporaneous")
  expect_error(MCMmodel(data, kernel = "dyn"), "partial matching")
  expect_error(MCMmodel(data, kernel = c("dynamic", "contemporaneous")),
               "exactly one")
  expect_error(MCMmodel(data, kernel = 1), "exactly one")

  dynamic_model <- MCMmodel(data, kernel = "dynamic")
  expect_identical(dynamic_model$meta_data$kernel, "dynamic")
  expect_named(dynamic_model$num_matrices,
               c("B", "Tau", "Kappa", "L_G", "D2"))
  expect_error(MCMmodel(data, n_latent = 1, kernel = "dynamic"),
               "observed-state")
  expect_output(print(dynamic_model), "Kernel: Stationary Dynamic MCMSEM")
  expect_output(summary(dynamic_model), "Nominal degrees of freedom: 1")
})

test_that("old-style model metadata defaults to the contemporaneous kernel", {
  data <- make_moment_summary()
  expect_false("kernel" %in% names(data$meta_data))

  model <- MCMmodel(data)
  model$meta_data$kernel <- NULL
  expect_identical(MCMSEM:::.model_kernel(model), "contemporaneous")
  expect_identical(MCMSEM:::.kernel_label(MCMSEM:::.model_kernel(model)),
                   "Contemporaneous Structural MCMSEM")
  expect_s4_class(model$copy(), "mcmmodelclass")

  fit <- MCMfit(
    model, data, compute_se = FALSE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.05
  )
  expect_identical(MCMSEM:::.result_kernel(fit), "contemporaneous")
  expect_identical(fit$model$meta_data$kernel, "contemporaneous")
})

test_that("dynamic matrices, starts, and bounds are editable", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMedit(model, "B", c(1, 2), "beta_Y_to_X")
  expect_true("beta_Y_to_X" %in% model$param_names)
  expect_equal(model$bounds["L", "beta_Y_to_X"], -0.98)
  model <- MCMedit(model, "lbound", "beta_Y_to_X", -0.5)
  model <- MCMedit(model, "start", "beta_Y_to_X", 0.15)
  expect_equal(model$bounds["L", "beta_Y_to_X"], -0.5)
  expect_equal(model$num_matrices$B[1, 2], 0.15)

  model <- MCMedit(model, "lbound", "B", -0.4)
  b_parameters <- vapply(model$param_coords, `[[`, character(1), 1) == "B"
  expect_true(all(model$bounds["L", b_parameters] == -0.4))
})

test_that("contemporaneous implied moments retain their numeric baseline", {
  data <- make_moment_summary()
  model <- suppressWarnings(MCMmodel(
    data, n_latent = 0, constrained_a = FALSE,
    kernel = "contemporaneous"
  ))
  values <- c(
    b2_1 = 0.12, b1_2 = -0.18, s1 = 1.3, s2 = 0.8,
    sk1 = 0.7, sk2 = -0.4, k1 = 4.2, k2 = 5.1
  )
  model <- set_dynamic_parameters(model, values)
  observed <- list(
    M2 = torch_tensor(data$M2), M3 = torch_tensor(data$M3),
    M4 = torch_tensor(data$M4)
  )
  matrices <- MCMSEM:::.get_torch_matrices(
    model, torch_device("cpu"), observed$M2, observed$M3, observed$M4,
    torch_float32()
  )
  slownecker <- jit_compile(MCMSEM:::.jit_funcs[["slownecker"]])
  predicted <- MCMSEM:::.get_predicted_matrices(
    matrices$.par_list, matrices$torch_masks, matrices$torch_maps,
    matrices$base_matrices, TRUE, TRUE, TRUE, FALSE, slownecker
  )
  expect_equal(as.matrix(predicted$M2), matrix(c(
    1.256646514, -0.132226139, -0.132226139, 0.806886017
  ), 2), tolerance = 2e-6)
  expect_equal(as.matrix(predicted$M3), matrix(c(
    0.655882955, -0.123577930, -0.123577930, -0.023747670,
    -0.123577930, -0.023747670, -0.023747670, -0.378989577
  ), 2, byrow = TRUE), tolerance = 2e-6)
  expect_identical(MCMSEM:::.model_kernel(model), "contemporaneous")
})
