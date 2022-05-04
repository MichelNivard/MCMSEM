test_that("MCMSEM works", {
  expect_type(MCMSEM(simulate_data(), bootstrap_iter=10), "double")
})
