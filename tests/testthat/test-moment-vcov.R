test_that("moment preparation records mean-corrected central-moment covariance", {
  set.seed(20260729)
  x <- as.data.frame(matrix(rnorm(2400), ncol = 2))
  names(x) <- c("X", "Y")
  data <- MCMdatasummary(
    x, scale_data = FALSE, prep_asymptotic_se = TRUE,
    use_skewness = TRUE, use_kurtosis = TRUE
  )
  Omega <- MCMSEM:::.dynamic_tensor_to_matrix(data$SE$S.m)
  expect_true(data$SE$computed)
  expect_identical(data$SE$representation, "raw_central_moments")
  expect_true(data$SE$influence_function_corrected)
  expect_identical(dim(Omega), c(12L, 12L))
  expect_equal(Omega, t(Omega), tolerance = 1e-7)
  expect_gte(min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values),
             -1e-7)

  path <- tempfile(fileext = ".mcmdata")
  on.exit(unlink(path), add = TRUE)
  MCMsavesummary(data, path)
  restored <- MCMdatasummary(path = path)
  expect_identical(restored$SE$representation, "raw_central_moments")
  expect_true(restored$SE$influence_function_corrected)
  expect_identical(restored$SE$covariance_scale,
                   "Var(sample moment vector)")
})
