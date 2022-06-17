test_that("MCMSEM works", {
  expect_type({
    data <- simulate_data()
    model <- MCMmodel(data)
    res <- MCMfit(data, compute_se=FALSE)
  }, "double")
})