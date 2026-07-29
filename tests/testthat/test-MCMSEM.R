test_that("a contemporaneous model still fits", {
  data <- make_moment_summary()
  model <- MCMmodel(data)
  result <- MCMfit(
    model, data, compute_se = FALSE,
    optimizers = "rprop", optim_iters = 1, learning_rate = 0.01
  )
  expect_s4_class(result, "mcmresultclass")
  expect_type(result$loss, "double")
  expect_identical(result$kernel, "contemporaneous")
})
