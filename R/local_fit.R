# Function to quadratically scale loss if estimates are out of bounds
.loss_power_bounds <- function(par, torch_bounds) {
  # Does not work for some reason, loss does not change, then jumps to NA
  lbound_check <- torch_less_equal(torch_cat(par), torch_bounds[["L"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  ubound_check <- torch_greater_equal(torch_cat(par), torch_bounds[["U"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  lbound_dist <- torch_square(torch_cat(par) - torch_bounds[["L"]])  # Distance from lbound
  ubound_dist <- torch_square(torch_cat(par) - torch_bounds[["U"]])  # Distance from ubound
  ldist_sum <- torch_sqrt(torch_sum(lbound_dist * lbound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  udist_sum <- torch_sqrt(torch_sum(ubound_dist * ubound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  # 1e-15 is added to the lines above to prevent NAN gradients due to sqrt(0) issue
  pow <- torch_sum(ldist_sum) + torch_sum(udist_sum)
  return(pow)
}

# Objective function
.torch_objective <- function(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis) {
  A <- torch_add(base_matrices[['A']], torch_sum(torch_mul(torch_maps[['A']], .par_list[['A']]), dim=3))
  Fm <- torch_add(base_matrices[['Fm']], torch_sum(torch_mul(torch_maps[['Fm']], .par_list[['Fm']]), dim=3))
  S <- torch_add(base_matrices[['S']], torch_sum(torch_mul(torch_maps[['S']], .par_list[['S']]), dim=3))

  if (use_skewness) {Sk <- torch_add(base_matrices[['Sk']], torch_sum(torch_mul(torch_maps[['Sk']], .par_list[['Sk']]), dim=3))}
  if (use_kurtosis) {K <- torch_add(torch_add(torch_mul(torch_matmul(torch_matmul(torch_sqrt(S), base_matrices[['K']]), .torch_kron(torch_sqrt(S), .torch_kron(torch_sqrt(S), torch_sqrt(S)))), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3)), base_matrices[['K2']])}

  # Rstyle: M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), S), torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), torch_transpose(Fm, 1, 2))
  if (use_skewness) {
    # Rstyle: M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
    M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), Sk), .torch_kron(torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2))), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
  }
  if (use_kurtosis) {
    # Rstyle: M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))
    M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), K), torch_matmul(.torch_kron(.torch_kron( torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))
  }
  ### Loss function
  # Rstyle: value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  if (use_skewness & use_kurtosis) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(M3, m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']])) + lossfunc(torch_mul(M4, m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  } else if (use_skewness) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(M3, m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']]))
  } else if (use_kurtosis) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(M4, m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  }
  if (use_bounds) {
    pow <- .loss_power_bounds(.par_list, torch_bounds)
    loss <- value * (2^pow)
    return(loss)
  } else {
    return(value)
  }
}

# Fit wrapper function
.torch_fit <- function(M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, lossfunc, return_history=FALSE, low_memory) {
  loss_hist <- NULL
  optim <- optim_rprop(.par_list,lr = learning_rate)
  if (low_memory) {gc(verbose=FALSE, full=TRUE)}
  for (i in 1:optim_iters[1]) {
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis)
    if (low_memory) {gc(verbose=FALSE, full=TRUE)}
    loss$backward()
    loss_hist <- c(loss_hist, loss$detach())
    optim$step()
    if (!(silent)) cat(paste0("\rloss ", as.numeric(loss), "           "))
  }
  calc_loss_torchfit <- function() {
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis)
    if (low_memory) {gc(verbose=FALSE, full=TRUE)}
    loss$backward()
    loss_hist <<- c(loss_hist, loss$detach())
    if (!(silent)) {cat(paste0("\rloss ", as.numeric(loss), "          "))}
    return(loss)
  }
  # Use lbfgs to get really close....
  learning_rate <- 1
  optim <- optim_lbfgs(.par_list,lr= learning_rate)
  for (i in 1:optim_iters[2]) {
    optim$step(calc_loss_torchfit)
  }
  if (!(silent)) {cat("\n")}
  if (return_history) {
    return(list(par=.par_list, loss_hist=loss_hist))
  } else {
    return(.par_list)
  }
}
