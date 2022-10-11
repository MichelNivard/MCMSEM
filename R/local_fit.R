# Function to quadratically scale loss if estimates are out of bounds
.loss_power_bounds <- function(par, torch_bounds, outofbounds_penalty) {
  # Does not work for some reason, loss does not change, then jumps to NA
  lbound_check <- torch_less_equal(torch_cat(par), torch_bounds[["L"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  ubound_check <- torch_greater_equal(torch_cat(par), torch_bounds[["U"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  lbound_dist <- torch_square(torch_cat(par) - torch_bounds[["L"]])  # Distance from lbound
  ubound_dist <- torch_square(torch_cat(par) - torch_bounds[["U"]])  # Distance from ubound
  ldist_sum <- torch_sqrt(torch_sum(lbound_dist * lbound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  udist_sum <- torch_sqrt(torch_sum(ubound_dist * ubound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  # 1e-15 is added to the lines above to prevent NAN gradients due to sqrt(0) issue
  pow <- torch_sum(ldist_sum) + torch_sum(udist_sum) * outofbounds_penalty
  return(pow)
}

# Calculate predicted M2, M3, M4 (depending on use_skewness, use_kurtosis) matrices
.get_predicted_matrices <- function(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, diag_s, low_memory) {
  # dev note: All the lines for computing these matrices (especially K, M2, M3 and M4) are unreadable due to their excessive length, I know.
  #           This is because it saves VRAM, any intermediate objects that are created here are also stored on the GPU and waste highly valuable space.
  #           I know this makes it a pain to modify antyhing, and I'm sorry for that, but it is something we have to deal with for now :(
  #           Then why not make the whole thing one line you ask? Well I tried that, and oddly enough that actually takes up even more memory...
  #           Don't ask me why but the current balance between storing intermediates and one-lining seems ideal.
  A <- torch_add(base_matrices[['A']], torch_sum(torch_mul(torch_maps[['A']], .par_list[['A']]), dim=3))
  Fm <- torch_add(base_matrices[['Fm']], torch_sum(torch_mul(torch_maps[['Fm']], .par_list[['Fm']]), dim=3))
  S <- torch_add(base_matrices[['S']], torch_sum(torch_mul(torch_maps[['S']], .par_list[['S']]), dim=3))

  if (use_skewness) {Sk <- torch_add(base_matrices[['Sk']], torch_sum(torch_mul(torch_maps[['Sk']], .par_list[['Sk']]), dim=3))}
  if (use_kurtosis) {
    # Rstyle:
    # (sqrt(S) %*% base_matrices[['K']] %*%  (sqrt(S) %o% torch_sqrt(S) %o% sqrt(S))) * base_matrices[['K2']] * torch_masks[['K']] + sum(torch_maps[['K']] * .par_list[['K']], dim=3)
    # Note this is purely for readability as this R code will not actually work with base-R matrices as it uses 3D tensors in sum(..., dim=3)
    ## This is still in testing phase:
    ## sqrts <- torch_sqrt(S)
    if (diag_s) {
      sqrts <- torch_sign(S) * torch_sqrt(torch_abs(S))
      skron <- .torch_kron(torch_diag(sqrts), .torch_kron(torch_diag(sqrts), torch_diag(sqrts)))
      K <- torch_add(torch_mul(torch_mul(torch_mul(torch_matmul(sqrts, base_matrices[['K']]), skron), base_matrices[['K2']]), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3))
    } else {
      # This is the 'true' version of K but is significantly more memory intensive
      sqrts <- torch_sign(S) * torch_sqrt(torch_abs(S))
      K <- torch_add(torch_mul(torch_mul(torch_matmul(torch_matmul(sqrts, base_matrices[['K']]), .torch_kron(sqrts, .torch_kron(sqrts, sqrts))), base_matrices[['K2']]), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3))
    }
  }

  # Rstyle: M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), S), torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), torch_transpose(Fm, 1, 2))
  if (use_skewness) {
    # Rstyle: M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
    M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), Sk), .torch_kron(torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2))), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
  }
  if (use_kurtosis) {
    # Rstyle: M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))
    if (low_memory) {
      # M4 <- .jit_slownecker$fn(.jit_slownecker$fn(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), K), torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), torch_transpose(Fm, 1, 2))
      fmkronrow <- torch_tensor(.torch_kron(.torch_kron(torch_sum(Fm, dim=1), torch_sum(Fm, dim=1)), torch_sum(Fm, dim=1)), device=Fm$device, dtype=torch_bool())
      M4 <- .jit_slownecker$fn(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), K), torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2), fmkronrow)
    } else {
      M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), K), torch_matmul(.torch_kron(.torch_kron( torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))
    }
  }
  pred_matrices <- if (use_kurtosis & use_skewness) {return(list(M2=M2, M3=M3, M4=M4))} else if (use_skewness) {return(list(M2=M2, M3=M3))} else if (use_kurtosis) {return(list(M2=M2, M4=M4))} else {return(list(M2=M2))}
}

.calc_loss <- function(lossfunc, pred_matrices, m2v_masks, M2.obs, M3.obs, M4.obs, use_skewness, use_kurtosis) {
  ### Loss function
  # Rstyle (in case of mse): value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  if (use_skewness & use_kurtosis) {
    value <- lossfunc(torch_mul(pred_matrices[['M2']], m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(pred_matrices[['M3']], m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']])) + lossfunc(torch_mul(pred_matrices[['M4']], m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  } else if (use_skewness) {
    value <- lossfunc(torch_mul(pred_matrices[['M2']], m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(pred_matrices[['M3']], m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']]))
  } else if (use_kurtosis) {
    value <- lossfunc(torch_mul(pred_matrices[['M2']], m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) + lossfunc(torch_mul(pred_matrices[['M4']], m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  } else {
    value <- lossfunc(torch_mul(pred_matrices[['M2']], m2v_masks[['m2']]))
  }
  return(value)
}

# Objective function
.torch_objective <- function(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis, outofbounds_penalty, diag_s, low_memory) {
  pred_matrices <- .get_predicted_matrices(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, diag_s, low_memory)
  value <- .calc_loss(lossfunc, pred_matrices, m2v_masks, M2.obs, M3.obs, M4.obs, use_skewness, use_kurtosis)
  if (use_bounds) {
    pow <- .loss_power_bounds(.par_list, torch_bounds, outofbounds_penalty)
    loss <- value * (2^pow)
    return(loss)
  } else {
    return(value)
  }
}

# Fit wrapper function
.torch_fit <- function(optimfunc, M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, lossfunc, return_history=FALSE, low_memory, outofbounds_penalty, diag_s,debug=FALSE, monitor_grads=FALSE) {
  loss_hist <- NULL
  optim <- optimfunc(.par_list,lr = learning_rate[1])
  #if (low_memory) {gc(verbose=FALSE, full=TRUE)}  # Superceeded by .jit_slowneckerproduct, turn this back on in case of memory issues
  for (i in seq_len(optim_iters[1])) {
    if (debug) {cat(paste0("Optimizer1:",i,"\n"))}
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis, outofbounds_penalty, diag_s, low_memory)
    #if (low_memory) {gc(verbose=FALSE, full=TRUE)}  # Superceeded by .jit_slowneckerproduct, turn this back on in case of memory issues
    loss$backward()
    if (monitor_grads) {
      for (matname in names(.par_list)) {if (.par_list[[matname]]$requires_grad) {if (any(as.logical(torch_isnan(torch_tensor(A, device=torch_device('cpu')))))) {
            warning(paste0("NaN gradients found in ", matname, ", MCMfit stopped early"))
              if (return_history) {pred_matrices <- .get_predicted_matrices(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, diag_s, low_memory);return(list(par=.par_list, loss_hist=loss_hist, pred_matrices=pred_matrices))
              } else {return(.par_list)}}}}}
    loss_hist <- c(loss_hist, loss$detach())
    optim$step()
    if (!(silent)) cat(paste0("\rloss ", as.numeric(loss), "           "))
  }
  calc_loss_torchfit <- function() {
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis, outofbounds_penalty, diag_s, low_memory)
    #if (low_memory) {gc(verbose=FALSE, full=TRUE)} # Superceeded by .jit_slowneckerproduct, turn this back on in case of memory issues
    loss$backward()
    loss_hist <<- c(loss_hist, loss$detach())
    if (!(silent)) {cat(paste0("\rloss ", as.numeric(loss), "          "))}
    return(loss)
  }
  # Use lbfgs to get really close....
  optim <- optim_lbfgs(.par_list,lr= learning_rate[2])
  for (i in seq_len(optim_iters[2])) {
    if (debug) {cat(paste0("Optimizer2:",i,"\n"))}
    optim$step(calc_loss_torchfit)
  }
  if (!(silent)) {cat("\n")}
  if (return_history) {
    pred_matrices <- .get_predicted_matrices(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, diag_s, low_memory)
    return(list(par=.par_list, loss_hist=loss_hist, pred_matrices=pred_matrices))
  } else {
    return(.par_list)
  }
}
