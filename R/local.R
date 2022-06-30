.m3m2v <- utils::getFromNamespace("M3.mat2vec","PerformanceAnalytics")
.m4m2v <- utils::getFromNamespace("M4.mat2vec","PerformanceAnalytics")
.m2m2v <- function(x){c(x[lower.tri(x,diag=T)])}
# Just in case M3.mat2vec and M4.mat2vec from Performanceanalytics ever changes I am leaving R alternatives for these functions here as translated from their C++ code
# Uncomment these, and expand them to make them readable when that happens
# .m3m2v <- function(x) {p <- nrow(x); M3vec <- rep(0, p * (p + 1) * (p + 2) / 6); iter <- 1; for (i in 0:(p-1)) {for (j in i:(p-1)) {for (k in j:(p-1)) {M3vec[iter] <- x[((i * p + j) * p + k)+1]; iter <- iter + 1}}}; return(M3vec)}
# .m4m2v <- function(x) {p <- nrow(x); M4vec <- rep(0, p * (p + 1) * (p + 2) * (p + 3) / 24); iter <- 1; for (i in 0:(p-1)) {for (j in i:(p-1)) {for (k in j:(p-1)) {for (l in k:(p-1)) {M4vec[iter] <- x[((i * p * p + j * p + k) * p + l) + 1]; iter <- iter + 1}}}}; return(M4vec)}

.gen_matrices <- function(par, n_p, n_f, base_value=0) {
  matrices <- list()
  ############## A
  A <- matrix(rep(base_value, (n_p+n_f)^2), n_p+n_f)
  for (i in 1:length(par[['a']])) {
    a_val <- par[['a']][i]
    A[((i-1) %% n_p)+n_f+1, floor((i-1)/n_p)+1] <- a_val
  }
  B <- matrix(rep(base_value, (n_p)^2), n_p)
  iter <- 1
  for (i in 1:n_p) {
    for (j in 1:n_p) {
      if (i > j) {
        B[i, j] <- par[['b']][iter]
        iter <- iter + 1
      } else if (i < j) {
        B[i, j] <- par[['b']][iter]
        iter <- iter + 1
      }
    }
  }
  A[(n_f+1):nrow(A), (n_f+1):ncol(A)] <- B
  matrices[['A']] <- A
  ############## Fm
  Fm <- matrix(rep(base_value, n_p * (n_p+n_f)), n_p)
  for (i in 1:n_p) {
    Fm[i, i+n_f] <- if (is.character(base_value)) "1" else 1
  }
  matrices[['Fm']] <- Fm
  ############## S
  S <- matrix(rep(base_value, (n_p+n_f)^2), n_p+n_f)
  for (i in 1:n_f) {
    S[i, i] <- if (is.character(base_value)) "1" else 1
  }
  for (i in 1:n_p) {
    S[i+n_f, i+n_f] <- par[['s']][i]
  }
  matrices[['S']] <- S
  ############## Sk
  Sk <- matrix(rep(base_value, (n_p+n_f)^3), n_p + n_f)
  for (i in 1:n_p) {
    Sk[i+n_f, i+n_f+(i-1+n_f)*(n_p+n_f)]  <- par[['sk']][i]
  }
  matrices[['Sk']] <- Sk
  ############## K
  K <- matrix(rep(base_value, (n_p+n_f)^4), n_p + n_f)
  if (!(is.character(base_value))) {
    K1_ref <- rnorm(2000)
    for (i in 1:(n_p + n_f - 1)) {
      K1_ref <- cbind(K1_ref, rnorm(2000))
    }
    K <- sqrt(S) %*% K %*% (sqrt(S) %x% sqrt(S) %x% sqrt(S))
    K1 <- K[ round(M4.MM(K1_ref)) != 0] <- 1
    matrices[['K1_ref']] <- round(M4.MM(K1_ref)) != 0
  }
  for (i in 1:n_p) {
    K[i+n_f, i+n_f + (i-1+n_f)*(n_p+n_f) + (i-1+n_f)*((n_p+n_f)^2)]  <- par[['k']][i]
  }
  matrices[['K']] <- K
  return(matrices)
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
  } else {
    stop("Oops, something went wrong. 3")
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
  if (use_skewness) {
    Sk <- torch_add(base_matrices[['Sk']], torch_sum(torch_mul(torch_maps[['Sk']], .par_list[['Sk']]), dim=3))
  }
  if (use_kurtosis) {
    K <- torch_add(torch_add(torch_mul(torch_matmul(torch_matmul(torch_sqrt(S), base_matrices[['K']]), .torch_kron(torch_sqrt(S), .torch_kron(torch_sqrt(S), torch_sqrt(S)))), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3)), base_matrices[['K2']])
  }

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
  } else {
    stop("Oops, something went wrong. 4")
  }
}

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

### pull it together to make std errors:
.std.err <- function(data, .par_list, use_skewness, use_kurtosis, torch_masks, torch_maps, base_matrices, m2v_masks, device) {
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
  } else {
    stop("Oops, something went wrong. 2")
  }
  Rm2vmasks <- list(
    m2=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m2']], device=torch_device("cpu"))))),
    m3=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m3']], device=torch_device("cpu"))))),
    m4=which(as.logical(as.numeric(torch_tensor(m2v_masks[['m4']], device=torch_device("cpu")))))
  )
  S.m <- apply(scale(data, center = T, scale = F),1, .t4crossprod, idx=idx, use_skewness=use_skewness, use_kurtosis=use_kurtosis) # 00:24

  S.m <- cov(t(S.m))/(n-1) # 1:01

  # weights matrix is based on diagonal, may be better behaved?
  W <- solve(diag(nrow(S.m)) * S.m) # 00:00.1
  par_vec <- c()
  par_to_list_coords <- list()
  coord_start <- 1
  for (i in names(.par_list)) {
    if (.par_list[[i]]$requires_grad) {
      current_vec <- as.numeric(torch_tensor(.par_list[[i]], device=torch_device("cpu")))
      par_vec <- c(par_vec, current_vec)
      par_to_list_coords[[i]] <- coord_start:(coord_start+(length(current_vec)-1))
      coord_start <- coord_start + length(current_vec)
    }
  }
  G <- jacobian(func = .jac.fn_torch,x = par_vec, .par_list=.par_list, par_to_list_coords=par_to_list_coords, torch_masks=torch_masks,
                torch_maps=torch_maps, base_matrices=base_matrices, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                Rm2vmasks=Rm2vmasks, device=device)

  Asycov <- solve(t(G)%*%W%*%G) %*% t(G)%*%W%*%S.m %*%W%*%G %*% solve(t(G)%*%W%*%G)
  se <- sqrt(2)*sqrt(diag(Asycov))
  return(se)

}

# wrapper function to make the code more R-like
.torch_kron <- function(a, b) {
  return(a$kron(b))
}

# Turn 1D index into 2D index
.r_1to2d_idx <- function(x, nrows) {
  row <- (x-1) %% nrows + 1
  col <- floor((x-1) / nrows) +1
  return(c(row, col))
}

.get_torch_coords <- function(model, device) {
  # parameter coordinates need to be translated as torch does not accept single intiger index for 2d matrix
  torch_coords <- list()
  for (i in 1:length(model$param_coords)) {
    i_par <- model$param_coords[[i]]
    i_mat_name <- i_par[[1]]
    i_coords <- i_par[[2]]
    i_mult <- i_par[[3]]
    i_mat <- model$num_matrices[[i_mat_name]]
    for (j in 1:length(i_coords)) {
      i_row_coords <- .r_1to2d_idx(i_coords[j], nrow(i_mat))[1]
      i_col_coords <- .r_1to2d_idx(i_coords[j], nrow(i_mat))[2]
      torch_coords <- append(torch_coords, list(list(
        mat_name=i_mat_name, row=i_row_coords, col=i_col_coords,
        par=i, mult=torch_tensor(i_mult[j], device=device))))
    }
  }
  return(torch_coords)
}

.torch_m2m2v <- function(x) {
  return(x[lower.tri(x, diag=TRUE)])
}
.torch_m3m2v <- function(x) {
  p <- x$shape[1]
  M3vec <- torch_empty(p * (p + 1) * (p + 2) / 6)
  iter <- 1
  for (i in 0:(p-1)) {
    for (j in i:(p-1)) {
      for (k in j:(p-1)) {
        coords <- .r_1to2d_idx(((i * p + j) * p + k)+1, p)
        M3vec[iter] <- x[coords[1], coords[2]]
        iter <- iter + 1
      }
    }
  }
  return(M3vec)
}

.torch_m4m2v <- function(x) {
  p <- x$shape[1]
  M4vec <- torch_empty(p * (p + 1) * (p + 2) * (p + 3) / 24)
  iter <- 1
  for (i in 0:(p-1)) {
    for (j in i:(p-1)) {
      for (k in j:(p-1)) {
        for (l in k:(p-1)) {
          coords <- .r_1to2d_idx(((i * p * p + j * p + k) * p + l) + 1, p)
          M4vec[iter] <- x[coords[1], coords[2]]
          iter <- iter + 1
        }
      }
    }
  }
  return(M4vec)
}


.torch_m2m2v_mask <- function(x, device, dtype) {
  mask <- torch_zeros_like(x, device=device, dtype=dtype)
  mask[lower.tri(x, diag=TRUE)] <- 1
  return(mask)
}
.torch_m3m2v_mask <- function(x, device, dtype) {
  mask <- torch_zeros_like(x, device=device, dtype=dtype)
  p <- mask$shape[1]
  for (i in 0:(p-1)) {
    for (j in i:(p-1)) {
      for (k in j:(p-1)) {
        coords <- .r_1to2d_idx(((i * p + j) * p + k)+1, p)
        mask[coords[1], coords[2]] <- 1
      }
    }
  }
  return(mask)
}

.torch_m4m2v_mask <- function(x, device, dtype) {
  mask <- torch_zeros_like(x, device=device, dtype=dtype)
  p <- mask$shape[1]
  for (i in 0:(p-1)) {
    for (j in i:(p-1)) {
      for (k in j:(p-1)) {
        for (l in k:(p-1)) {
          coords <- .r_1to2d_idx(((i * p * p + j * p + k) * p + l) + 1, p)
          mask[coords[1], coords[2]] <- 1
        }
      }
    }
  }
  return(mask)
}

.get_torch_matrices <- function(model, device, M2.obs, M3.obs, M4.obs, torch_dtype) {
    torch_coords <- .get_torch_coords(model, device)

  torch_matrices <- list(
    A=torch_tensor(model$num_matrices[["A"]], device=device, dtype=torch_dtype),
    Fm= torch_tensor(model$num_matrices[["Fm"]], device=device, dtype=torch_dtype),
    S=torch_tensor(model$num_matrices[["S"]], device=device, dtype=torch_dtype),
    Sk=torch_tensor(model$num_matrices[["Sk"]], device=device, dtype=torch_dtype),
    K=torch_tensor(model$num_matrices[["K1_ref"]]+1-1, device=device, dtype=torch_dtype),
    diag_n_p=torch_tensor(torch_diagflat(rep(1, model$meta_data$n_phenotypes + model$meta_data$n_confounding)), device=device, dtype=torch_dtype)
  )
  param_list <- list(A=c(), Fm= c(), S=c(), Sk=c(), K=c())
  torch_maps <- list(A=list(), Fm= list(), S=list(), Sk=list(), K=list())
  torch_masks <- list(
    A=torch_ones_like(torch_matrices[["A"]], dtype=torch_dtype, device=device),
    Fm= torch_ones_like(torch_matrices[["Fm"]], dtype=torch_dtype, device=device),
    S=torch_ones_like(torch_matrices[["S"]], dtype=torch_dtype, device=device),
    Sk=torch_ones_like(torch_matrices[["Sk"]], dtype=torch_dtype, device=device),
    K=torch_ones_like(torch_matrices[["K"]], dtype=torch_dtype, device=device)
  )
  for (i in torch_coords) {
    new_mat <- torch_zeros_like(torch_matrices[[i$mat_name]], device=device, dtype=torch_dtype)
    new_mat[i$row, i$col] <- 1 * i$mult
    torch_masks[[i$mat_name]][i$row, i$col] <- 0.0
    param_name <- model$named_matrices[[i$mat_name]][i$row, i$col]
    if (param_name %in% names(torch_maps[[i$mat_name]])) {
      torch_maps[[i$mat_name]][[param_name]] <- torch_maps[[i$mat_name]][[param_name]]  + new_mat
    } else {
      torch_maps[[i$mat_name]][[param_name]] <- new_mat
      param_list[[i$mat_name]] <- c(param_list[[i$mat_name]], model$start_values[param_name])
    }
  }

  n_p <- model$meta_data$n_phenotypes + model$meta_data$n_confounding
  for (i in 1:model$meta_data$n_confounding) {
    torch_masks[['K']][i, i + (i-1)*(n_p) + (i-1)*((n_p)^2)]  <- 0
  }
  # 3D tensors defining locations of paramters that require grad
  torch_maps <- list(
    A=if (length(torch_maps[['A']]) > 0) {torch_dstack(torch_maps[['A']])} else {torch_zeros_like(torch_matrices[["A"]], device=device, dtype=torch_dtype)},
    Fm=if (length(torch_maps[['Fm']]) > 0) {torch_dstack(torch_maps[['Fm']])} else {torch_zeros_like(torch_matrices[["Fm"]], device=device, dtype=torch_dtype)},
    S=if (length(torch_maps[['S']]) > 0) {torch_dstack(torch_maps[['S']])} else {torch_zeros_like(torch_matrices[["S"]], device=device, dtype=torch_dtype)},
    Sk=if (length(torch_maps[['Sk']]) > 0) {torch_dstack(torch_maps[['Sk']])} else {torch_zeros_like(torch_matrices[["Sk"]], device=device, dtype=torch_dtype)},
    K=if (length(torch_maps[['K']]) > 0) {torch_dstack(torch_maps[['K']])} else {torch_zeros_like(torch_matrices[["K"]], device=device, dtype=torch_dtype)}
  )
  # Reshape to 3D is necessary as we sum over 3rd axis later on
  for (i in names(torch_maps)) {
    if (length(torch_maps[[i]]$shape) == 2) {
      shape <- as.numeric(torch_maps[[i]]$shape)
      torch_maps[[i]] <- torch_reshape(torch_maps[[i]], c(shape[1], shape[2], 1))
    }
  }
  K2 <- torch_zeros_like(torch_matrices[['K']], device=device, dtype=torch_dtype)
  for (i in 1:model$meta_data$n_confounding) {
    K2[i, i + (i-1)*(n_p) + (i-1)*((n_p)^2)]  <- 3
  }
  base_matrices <- list(
    A=torch_mul(torch_tensor(model$num_matrices[["A"]], device=device, dtype=torch_dtype), torch_masks[['A']]),
    Fm=torch_mul(torch_tensor(model$num_matrices[["Fm"]], device=device, dtype=torch_dtype), torch_masks[['Fm']]),
    S=torch_mul(torch_tensor(model$num_matrices[["S"]], device=device, dtype=torch_dtype), torch_masks[['S']]),
    Sk=torch_mul(torch_tensor(model$num_matrices[["Sk"]], device=device, dtype=torch_dtype), torch_masks[['Sk']]),
    K=torch_tensor(model$num_matrices[["K1_ref"]]+1-1, device=device, dtype=torch_dtype),
    K2=K2,
    diag_n_p=torch_tensor(torch_diagflat(rep(1, model$meta_data$n_phenotypes + model$meta_data$n_confounding)), device=device, dtype=torch_dtype)
  )
  # Zeros can produce NAN gradients, therefore set values in S matrices to very low values
  base_matrices[['S']] <- base_matrices[['S']] + torch_tensor(1e-16, device=device)
  # If the matrix (A, Fm, S, Sk, or K) does not have free parameters, use a constant 1 to ensure code will always work
  .par_list <- list(
    A=if (length(param_list[['A']]) > 0) {torch_tensor(as.numeric(param_list[['A']]), requires_grad = TRUE, device=device, dtype=torch_dtype)} else {torch_tensor(1, device=device, dtype=torch_dtype)},
    Fm=if (length(param_list[['Fm']]) > 0) {torch_tensor(as.numeric(param_list[['Fm']]), requires_grad = TRUE, device=device, dtype=torch_dtype)} else {torch_tensor(1, device=device, dtype=torch_dtype)},
    S=if (length(param_list[['S']]) > 0) {torch_tensor(as.numeric(param_list[['S']]), requires_grad = TRUE, device=device, dtype=torch_dtype)} else {torch_tensor(1, device=device, dtype=torch_dtype)},
    Sk=if (length(param_list[['Sk']]) > 0) {torch_tensor(as.numeric(param_list[['Sk']]), requires_grad = TRUE, device=device, dtype=torch_dtype)} else {torch_tensor(1, device=device, dtype=torch_dtype)},
    K=if (length(param_list[['K']]) > 0) {torch_tensor(as.numeric(param_list[['K']]), requires_grad = TRUE, device=device, dtype=torch_dtype)} else {torch_tensor(1, device=device, dtype=torch_dtype)}
  )
  torch_bounds <- list(L=list(), U=list())
  for (i in names(.par_list)) {
    if (.par_list[[i]]$requires_grad) {
      torch_bounds[["L"]][[i]] <- torch_tensor(as.numeric(model$bounds["L", names(param_list[[i]])]), device=device, dtype=torch_dtype)
      torch_bounds[["U"]][[i]] <- torch_tensor(as.numeric(model$bounds["U", names(param_list[[i]])]), device=device, dtype=torch_dtype)
    } else {
      # Non-grad parameters also need bounds in order for the shapes to match
      torch_bounds[["L"]][[i]] <- torch_ones_like(.par_list[[i]], device=device, dtype=torch_dtype) - 100
      torch_bounds[["U"]][[i]] <- torch_tensor(.par_list[[i]], device=device, dtype=torch_dtype) + 100
    }
  }
  torch_bounds[['L']] <- torch_cat(torch_bounds[['L']])
  torch_bounds[['U']] <- torch_cat(torch_bounds[['U']])
  m2v_masks <- list(
    m2=.torch_m2m2v_mask(M2.obs, device=device, dtype=torch_dtype),
    m3=.torch_m3m2v_mask(M3.obs, device=device, dtype=torch_dtype),
    m4=.torch_m4m2v_mask(M4.obs, device=device, dtype=torch_dtype)
  )
  torch_masks <- list(K=torch_masks[['K']])
  return(list(
    m2v_masks=m2v_masks,
    param_list=param_list,
    torch_bounds=torch_bounds,
    torch_masks=torch_masks,
    torch_maps=torch_maps,
    base_matrices=base_matrices,
    .par_list=.par_list))
}
