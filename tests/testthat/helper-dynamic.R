make_moment_summary <- function(p = 2L, names = LETTERS[seq_len(p)], n = 10000L) {
  M2 <- diag(seq(1.5, 1.5 + 0.2 * (p - 1L), length.out = p), p)
  D3 <- array(0, rep(p, 3L))
  D4 <- array(0, rep(p, 4L))
  for (i in seq_len(p)) {
    D3[matrix(rep(i, 3L), nrow = 1)] <- (-1)^i * (0.5 + i / 10)
    D4[matrix(rep(i, 4L), nrow = 1)] <- 1 + i / 2
  }
  M4 <- MCMSEM:::.cumulant4_to_raw4(D4, M2)
  mcmdataclass(
    meta_data = list(
      scale_data = FALSE, data_was_scaled = TRUE, weighted = FALSE,
      ncol = p, colnames = names, N = n, weightsum = n
    ),
    M2 = M2,
    M3 = MCMSEM:::.dynamic_array_to_mcm(D3),
    M4 = MCMSEM:::.dynamic_array_to_mcm(M4),
    SE = list(computed = FALSE, S.m = NULL, idx = NULL)
  )
}

set_dynamic_parameters <- function(model, values) {
  for (name in names(values)) model <- MCMedit(model, "start", name, values[[name]])
  model
}

make_population_summary <- function(model, n = 10000L) {
  implied <- MCMSEM:::.dynamic_implied_moments_base(model)
  mcmdataclass(
    meta_data = list(
      scale_data = model$meta_data$scale_data,
      data_was_scaled = TRUE,
      weighted = FALSE,
      ncol = model$meta_data$n_phenotypes,
      colnames = model$meta_data$original_colnames,
      N = n,
      weightsum = n
    ),
    M2 = implied$M2, M3 = implied$M3, M4 = implied$M4,
    SE = list(computed = FALSE, S.m = NULL, idx = NULL)
  )
}

prepare_population_vcov <- function(data, Omega = NULL) {
  n_moments <- unname(MCMmomentcount(data)[["total"]])
  if (is.null(Omega)) Omega <- diag(n_moments) / data$meta_data$N
  data$SE <- list(
    computed = TRUE,
    S.m = torch_tensor(Omega),
    idx = list(idx = seq_len(n_moments)),
    representation = "raw_central_moments",
    influence_function_corrected = TRUE,
    covariance_scale = "Var(sample moment vector)"
  )
  data
}

simulate_dynamic_cross_section <- function(n, B, tau_kind = c("gamma", "chisq"),
                                           Psi_G, burnin = 150L, seed = 1L) {
  set.seed(seed)
  p <- nrow(B)
  state <- matrix(0, n, p)
  for (time in seq_len(burnin)) {
    innovations <- cbind(
      (rgamma(n, shape = 4) - 4) / 2,
      -(rchisq(n, df = 5) - 5) / sqrt(10)
    )
    state <- state %*% t(B) + innovations
  }
  gaussian <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Psi_G)
  out <- as.data.frame(state + gaussian)
  names(out) <- c("X", "Y")
  out
}
