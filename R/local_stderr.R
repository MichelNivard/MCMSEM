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

# Make pseudo obs for standard errors:
.t4crossprod <- function(x, idx, use_skewness, use_kurtosis){
  x2 <- x %o% x
  x3 <- x2 %x% x
  if (use_kurtosis) {
    x4 <- x3 %x% x
  }
  if (use_skewness & use_kurtosis) {
    return(c(as.vector(t(x2)), as.vector(t(x3)), as.vector(t(x4)))[idx])
  } else if (use_skewness) {
    return(c(as.vector(t(x2)), as.vector(t(x3)))[idx])
  } else if (use_kurtosis) {
    return(c(as.vector(t(x2)), as.vector(t(x4)))[idx])
  }
}

######## Compute Jacobian:
.jac.fn_torch <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device) {
  for (i in names(par_to_list_coords)) {
    .par_list[[i]] <- torch_tensor(par_vec[par_to_list_coords[[i]]], device=device)
  }
  A <- torch_add(base_matrices[['A']], torch_sum(torch_mul(torch_maps[['A']], .par_list[['A']]), dim=3))
  Fm <- torch_add(base_matrices[['Fm']], torch_sum(torch_mul(torch_maps[['Fm']], .par_list[['Fm']]), dim=3))
  S <- torch_add(base_matrices[['S']], torch_sum(torch_mul(torch_maps[['S']], .par_list[['S']]), dim=3))

  if (use_skewness) {Sk <- torch_add(base_matrices[['Sk']], torch_sum(torch_mul(torch_maps[['Sk']], .par_list[['Sk']]), dim=3))}
  if (use_kurtosis) {K <- torch_add(torch_add(torch_mul(torch_matmul(torch_matmul(torch_sqrt(S), base_matrices[['K']]), .torch_kron(torch_sqrt(S), .torch_kron(torch_sqrt(S), torch_sqrt(S)))), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3)), base_matrices[['K2']])}

  diag_n_p <- base_matrices[['diag_n_p']]
  ###### Compute the observed cov, cosk, and cokurt matrices #################
  diag_np_a_inv <- torch_inverse(diag_n_p - A)
  diag_np_a_inv_t <- torch_transpose(torch_inverse(diag_n_p - A), 1, 2)

  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), S), diag_np_a_inv_t), torch_transpose(Fm, 1, 2))
  M2 <- as.matrix(torch_tensor(M2, device=torch_device("cpu")))
  if (use_skewness) {
    M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), Sk), .torch_kron(diag_np_a_inv_t, diag_np_a_inv_t)), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
    M3 <- as.matrix(torch_tensor(M3, device=torch_device("cpu")))
  }
  if (use_kurtosis) {
    M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), K), torch_matmul(.torch_kron(.torch_kron(diag_np_a_inv_t, diag_np_a_inv_t), diag_np_a_inv_t), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))
    M4 <- as.matrix(torch_tensor(M4, device=torch_device("cpu")))
  }
  if (use_skewness & use_kurtosis) {
    return(c(M2[Rm2vmasks[['m2']]], M3[Rm2vmasks[['m3']]], M4[Rm2vmasks[['m4']]]))
  } else if (use_skewness) {
    return(c(M2[Rm2vmasks[['m2']]],M3[Rm2vmasks[['m3']]]))
  } else if (use_kurtosis) {
    return(c(M2[Rm2vmasks[['m2']]],M4[Rm2vmasks[['m4']]]))
  }
}

# Jacobian wrapper to force garbage collect after .jac.fn_torch_ call
.jac.fn_torch_lowmem <- function(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device) {
  val <- .jac.fn_torch(par_vec, .par_list, par_to_list_coords, torch_masks, torch_maps, base_matrices, use_skewness, use_kurtosis, Rm2vmasks, device)
  gc(verbose=FALSE, full=TRUE)
  return(val)
}

### pull it together to make std errors:
.std.err <- function(data, .par_list, use_skewness, use_kurtosis, torch_masks, torch_maps, base_matrices, m2v_masks, device, low_memory) {
  # .std.err <- function(data, .par_list, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_skewness, use_kurtosis){
  n <- nrow(data)
  # if there is too much data for spoeedly opperation, sample 100000 observations to base this on
  if(n > 100000){
    samp <- sample(1:n,100000,F)
    data <- data[samp, ]
  }

  # observed cov between pseudo obsertvations ovver n-1 gets us cov betwene moments moments
  dim2 <- data[1, ] %o% data[1, ]
  dim2locs <- .dimlocations(nrow(t(dim2)), dims=2)
  dim3 <- dim2 %x% data[1, ]
  dim3locs <- .dimlocations(nrow(t(dim3)), dims=3)
  dim4 <- dim3 %x% data[1, ]
  dim4locs <- .dimlocations(nrow(t(dim4)), dims=4)
  if (use_skewness & use_kurtosis) {
    idx <- c(dim2locs, length(as.vector(dim2)) + dim3locs,  length(as.vector(dim2)) + length(as.vector(dim3)) + dim4locs)
  } else if (use_skewness) {
    idx <- c(dim2locs, length(as.vector(dim2)) + dim3locs)
  } else if (use_kurtosis) {
    idx <- c(dim2locs, length(as.vector(dim2)) + dim4locs)
  }
  Rm2vmasks <- list(
    m2=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m2']], device=torch_device("cpu"))))),
    m3=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m3']], device=torch_device("cpu"))))),
    m4=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m4']], device=torch_device("cpu")))))
  )
  S.m <- apply(scale(data, center = T, scale = F),1, .t4crossprod, idx=idx, use_skewness=use_skewness, use_kurtosis=use_kurtosis) # 00:24

  # S.m <- cov(t(S.m))/(n-1)
  # Replace this with torch_cov(S.m) / (n - 1)  once that is implemented in torch for R
  S.m <- torch_transpose(torch_tensor(S.m, device=torch_device("cpu")), 1, 2)
  S.m <- torch_subtract(S.m, torch_reshape(torch_mean(S.m, dim=1), c(1, S.m$shape[2])))
  S.m <- torch_matmul(torch_transpose(S.m, 1, 2), S.m)  / (S.m$shape[1] - 1) / (n - 1)

  # weights matrix is based on diagonal, may be better behaved?
  S.m <- as.matrix(S.m)
  W <- solve(diag(nrow(S.m)) * S.m) # 00:00.1
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

  if (low_memory) {func <- .jac.fn_torch_lowmem} else {func <- .jac.fn_torch}
  G <- jacobian(func = func,x = par_vec, .par_list=.par_list, par_to_list_coords=par_to_list_coords, torch_masks=torch_masks,
                torch_maps=torch_maps, base_matrices=base_matrices, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                Rm2vmasks=Rm2vmasks, device=device)

  Asycov <- solve(t(G)%*%W%*%G) %*% t(G)%*%W%*%S.m %*%W%*%G %*% solve(t(G)%*%W%*%G)
  se <- sqrt(2)*sqrt(diag(Asycov))
  return(se)
}