.dimlocations <- function(p, dims=2) {
  if (dims == 2) vect <- rep(0, p * (p + 1)  / 3)
  if (dims == 3) vect <- rep(0, p * (p + 1) * (p + 2) / 6)
  if (dims == 4) vect <- rep(0, p * (p + 1) * (p + 2) * (p + 3) / 24)
  iter <- 1
  for (i in 0:(p-1)) {
    for (j in i:(p-1)) {
      if (dims == 2) {
        vect[iter] <- ((i * p + j))+1; iter <- iter+1
      } else {
        for (k in j:(p-1)) {
          if (dims == 3) {
            vect[iter] <- ((i * p + j) * p + k)+1; iter <- iter+1
          } else {
            for (l in k:(p-1)){
              vect[iter] <- (((i * p * p + j * p + k) * p + l) + 1); iter <- iter + 1
          }}}}}}
  return(vect)
}

######## Compute Jacobian:
.jac.fn_torch <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, m2vmasks1d, device, diag_s, low_memory, .jit_slownecker) {
  
 

  for (i in names(par_to_list_coords)) {
    .par_list[[i]] <- torch_tensor(par_vec[par_to_list_coords[[i]]], device=device)
  }
  pred_matrices <- .get_predicted_matrices(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis,diag_s, low_memory, .jit_slownecker)
  print(torch_flatten(pred_matrices[['M2']]))
  print(m2vmasks1d[['m2']])
  if (use_skewness & use_kurtosis) {
    return(as.numeric(torch_tensor(torch_hstack(list(torch_flatten(pred_matrices[['M2']])[m2vmasks1d[['m2']]], torch_flatten(pred_matrices[['M3']])[m2vmasks1d[['m3']]], torch_flatten(pred_matrices[['M4']])[m2vmasks1d[['m4']]])), device=torch_device('cpu'))))
  } else if (use_skewness) {
    return(as.numeric(torch_tensor(torch_hstack(list(torch_flatten(pred_matrices[['M2']])[m2vmasks1d[['m2']]], torch_flatten(pred_matrices[['M3']])[m2vmasks1d[['m3']]])), device=torch_device('cpu'))))
  } else if (use_kurtosis) {
    return(as.numeric(torch_tensor(torch_hstack(list(torch_flatten(pred_matrices[['M2']])[m2vmasks1d[['m2']]],  torch_flatten(pred_matrices[['M4']])[m2vmasks1d[['m4']]])), device=torch_device('cpu'))))
  } else {
    return(as.numeric(torch_tensor(torch_flatten(pred_matrices[['M2']])[m2vmasks1d[['m2']]], device=torch_device('cpu'))))
  }
}

# Jacobian wrapper to force garbage collect after .jac.fn_torch, superceeded by .jit_slowneckerproduct, re-enable this in case of memory issues
.jac.fn_torch_lowmem <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, m2vmasks1d, device, diag_s, low_memory, .jit_slownecker) {
  val <- .jac.fn_torch(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, m2vmasks1d, device, diag_s, low_memory, .jit_slownecker)
  cuda_empty_cache()
  return(val)
}

### pull it together to make std errors:
.std.err <- function(data, .par_list, use_skewness, use_kurtosis, torch_masks, torch_maps, base_matrices, m2v_masks, device, low_memory, diag_s, jacobian_method, debug) {
  # Temporary fix to check if forcing CPU leads to more consistency
  if ((device == torch_device('cpu')) & (.par_list[['A']]$is_cuda)) {
    if (debug) {cat(" - Reassigning tensors to CPU device\n")}
    # If CPU is used for SE but not for optimization, port all matrices to CPU first
    for (i in names(.par_list)) {.par_list[[i]] <- torch_tensor(.par_list[[i]], device=torch_device('cpu'), requires_grad = .par_list[[i]]$requires_grad)}
    for (i in names(torch_masks)) {torch_masks[[i]] <- torch_tensor(torch_masks[[i]], device=torch_device('cpu'))}
    for (i in names(torch_maps)) {torch_maps[[i]] <- torch_tensor(torch_maps[[i]], device=torch_device('cpu'))}
    for (i in names(base_matrices)) {base_matrices[[i]] <- torch_tensor(base_matrices[[i]], device=torch_device('cpu'))}
    for (i in names(m2v_masks)) {m2v_masks[[i]] <- torch_tensor(m2v_masks[[i]], device=torch_device('cpu'))}
  }

  # idx is pre-determined in MCMdatasummary() depending on the number of variables
  # Note there is probably a cleaner way to do this as these numbers only depend on number of columns..
  # If you want to be a contributor, here's your chance.

  par_vec <- NULL
  par_to_list_coords <- list()
  coord_start <- 1
  .par_list_grad_only <- list()
  if (debug) {cat(" - Reformatting parameters for jacobian\n")}
  for (i in names(.par_list)) {
    if (.par_list[[i]]$requires_grad) {
      current_vec <- as.numeric(torch_tensor(.par_list[[i]], device=torch_device("cpu")))
      par_vec <- c(par_vec, current_vec)
      par_to_list_coords[[i]] <- coord_start:(coord_start+(length(current_vec)-1))
      coord_start <- coord_start + length(current_vec)
      .par_list_grad_only[[i]] <- .par_list[[i]]
    } else {
      .par_list_grad_only[[i]] <- NULL
    }
  }
  if (debug) {cat(" - Reformatting masks for jacobian\n")}
  m2vmasks1d <- list(
    m2=torch_where(torch_flatten(torch_transpose(m2v_masks[['m2']], 1, 2))==1)[[1]] + torch_tensor(1, dtype=torch_long(), device=device),
    m3=torch_where(torch_flatten(torch_transpose(m2v_masks[['m3']], 1, 2))==1)[[1]] + torch_tensor(1, dtype=torch_long(), device=device),
    m4=torch_where(torch_flatten(torch_transpose(m2v_masks[['m4']], 1, 2))==1)[[1]] + torch_tensor(1, dtype=torch_long(), device=device)
  )
  # Slownecker function compiled when .std.err is called
  slowneckerfun <- if (low_memory > 2) {'slowernecker'}  else {'slownecker'}
  .jit_slownecker <- jit_compile(.jit_funcs[[slowneckerfun]])
  if (low_memory > 1) {func <- .jac.fn_torch_lowmem} else {func <- .jac.fn_torch}
  if (debug) {cat(" - Calculating jacobian\n")}
  G <- torch_tensor(jacobian(func = func,
                x = par_vec, method = jacobian_method, .par_list=.par_list, par_to_list_coords=par_to_list_coords, torch_masks=torch_masks,
                torch_maps=torch_maps, base_matrices=base_matrices, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                m2vmasks1d=m2vmasks1d, device=device, diag_s=diag_s, low_memory=low_memory, .jit_slownecker=.jit_slownecker),
                    device=torch_device('cpu'))
  # Note all these matrices used for Asycov are on CPU as this operation is only performe once, and the matrices involved are extremely large
  # weights matrix is based on diagonal, may be better behaved?
  # RStyle: W <- solve(diag(nrow(S.m)) * S.m)
  # Rstyle: Asycov <- solve(t(G)%*%W%*%G) %*% t(G)%*%W%*%S.m %*%W%*%G %*% solve(t(G)%*%W%*%G)
  if (debug) {cat(" - Calculating Asycov\n")}
  if (use_kurtosis & use_skewness) {
    # W is no longer assigned, and replaced by torch_inverse(torch_eye(nrow(S.m)) *S.m) to prevent unnecessary memory use
    # Additionally, data$SE$S.m is used to prevent creating a local copy of S.m
    Asycov <- torch_matmul(torch_matmul(torch_matmul(torch_inverse(torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(data$SE$S.m))*data$SE$S.m)), G)), torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(data$SE$S.m))*data$SE$S.m)), data$SE$S.m)), torch_matmul(torch_inverse(torch_eye(nrow(data$SE$S.m))*data$SE$S.m), G)), torch_inverse(torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(data$SE$S.m))*data$SE$S.m)), G)))
  } else {
    if (use_skewness) {
      idx <- data$SE$idx$idx_nokurt
    } else if (use_kurtosis) {
      idx <- data$SE$idx$idx_noskew
    } else {
      idx <- data$SE$idx$idx_nokurt_noskew
    }
    S.m <- data$SE$S.m[idx, idx]
    Asycov <- torch_matmul(torch_matmul(torch_matmul(torch_inverse(torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(S.m)) *S.m)), G)), torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(S.m)) *S.m)), S.m)), torch_matmul(torch_inverse(torch_eye(nrow(S.m)) *S.m), G)), torch_inverse(torch_matmul(torch_matmul(torch_transpose(G, 1, 2), torch_inverse(torch_eye(nrow(S.m)) *S.m)), G)))
  }
  se <- torch_sqrt(torch_diag(Asycov))
  return(as.numeric(se))
}
