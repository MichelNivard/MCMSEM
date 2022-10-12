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
.jac.fn_torch <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device, diag_s, low_memory, .jit_slownecker) {
  for (i in names(par_to_list_coords)) {
    .par_list[[i]] <- torch_tensor(par_vec[par_to_list_coords[[i]]], device=device)
  }
  pred_matrices <- .get_predicted_matrices(.par_list, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis,diag_s, low_memory, .jit_slownecker)
  pred_matrices[['M2']] <- as.matrix(torch_tensor(pred_matrices[['M2']], device=torch_device("cpu")))
  if (use_skewness) {
    pred_matrices[['M3']] <- as.matrix(torch_tensor(pred_matrices[['M3']], device=torch_device("cpu")))
  }
  if (use_kurtosis) {
    pred_matrices[['M4']] <- as.matrix(torch_tensor(pred_matrices[['M4']], device=torch_device("cpu")))
  }
  if (use_skewness & use_kurtosis) {
    return(c(pred_matrices[['M2']][Rm2vmasks[['m2']]], pred_matrices[['M3']][Rm2vmasks[['m3']]], pred_matrices[['M4']][Rm2vmasks[['m4']]]))
  } else if (use_skewness) {
    return(c(pred_matrices[['M2']][Rm2vmasks[['m2']]], pred_matrices[['M3']][Rm2vmasks[['m3']]]))
  } else if (use_kurtosis) {
    return(c(pred_matrices[['M2']][Rm2vmasks[['m2']]], pred_matrices[['M4']][Rm2vmasks[['m4']]]))
  } else {
    return(as.vector(pred_matrices[['M2']][Rm2vmasks[['m2']]]))
  }
}

# Jacobian wrapper to force garbage collect after .jac.fn_torch, superceeded by .jit_slowneckerproduct, re-enable this in case of memory issues
# .jac.fn_torch_lowmem <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device, diag_s, low_memory) {val <- .jac.fn_torch(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device, diag_s, low_memory); gc(verbose=FALSE, full=TRUE); return(val)}

### pull it together to make std errors:
.std.err <- function(data, .par_list, use_skewness, use_kurtosis, torch_masks, torch_maps, base_matrices, m2v_masks, device, low_memory, diag_s, jacobian_method) {
  # idx is pre-determined in MCMdatasummary() depending on the number of variables
  # Note there is probably a cleaner way to do this as these numbers only depend on number of columns..
  # If you want to be a contributor, here's your chance.
  if (use_kurtosis & use_skewness) {
    idx <- data$SE$idx$idx
  } else if (use_skewness) {
    idx <- data$SE$idx$idx_nokurt
  } else if (use_kurtosis) {
    idx <- data$SE$idx$idx_noskew
  } else {
    idx <- data$SE$idx$idx_nokurt_noskew
  }
  S.m <- data$SE$S.m[idx, idx]  # S.m is pre-computed in MCMdatasummary()
  # weights matrix is based on diagonal, may be better behaved?
  W <- solve(diag(nrow(S.m)) * S.m)
  par_vec <- NULL
  par_to_list_coords <- list()
  coord_start <- 1
  .par_list_grad_only <- list()
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

  Rm2vmasks <- list(
    m2=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m2']], device=torch_device("cpu"))))),
    m3=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m3']], device=torch_device("cpu"))))),
    m4=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m4']], device=torch_device("cpu")))))
  )
  # Slownecker function compiled when .std.err is called
  .jit_slownecker <- jit_compile(.jit_funcs[['slownecker']])
  # if (low_memory) {func <- .jac.fn_torch_lowmem} else {func <- .jac.fn_torch}  # superceeded by .jit_slowneckerproduct, turn this back on in case of memory issues
  G <- jacobian(func = .jac.fn_torch,  # If the previous line is re-enabled, replace with: G <- jacobian(func = func,
                x = par_vec, method = jacobian_method, .par_list=.par_list, par_to_list_coords=par_to_list_coords, torch_masks=torch_masks,
                torch_maps=torch_maps, base_matrices=base_matrices, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                Rm2vmasks=Rm2vmasks, device=device, diag_s=diag_s, low_memory=low_memory, .jit_slownecker=.jit_slownecker)
  list2env(list(G=G), .GlobalEnv) 
  list2env(list(W=W), .GlobalEnv) 
  Asycov <- solve(t(G)%*%W%*%G) %*% t(G)%*%W%*%S.m %*%W%*%G %*% solve(t(G)%*%W%*%G)
  se <- sqrt(2)*sqrt(diag(Asycov))
  return(se)
}
