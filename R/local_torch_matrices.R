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
  # Data types for torch matrices
  dtypes <- list(
    A=torch_dtype,
    Fm=torch_dtype,
    S=torch_dtype,
    Sk=torch_dtype,
    K=torch_dtype,
    diag_n_p=torch_dtype,
    bounds=torch_dtype
  )
  # parameter coordinates need to be translated as torch does not accept single intiger index for 2d matrix
  torch_coords <- list()
  for (i in seq_along(model$param_coords)) {
    i_par <- model$param_coords[[i]]
    i_mat_name <- i_par[[1]]
    i_coords <- i_par[[2]]
    i_mult <- i_par[[3]]
    i_mat <- model$num_matrices[[i_mat_name]]
    for (j in seq_along(i_coords)) {
      i_row_coords <- .r_1to2d_idx(i_coords[j], nrow(i_mat))[1]
      i_col_coords <- .r_1to2d_idx(i_coords[j], nrow(i_mat))[2]
      torch_coords <- append(torch_coords, list(list(
        mat_name=i_mat_name, row=i_row_coords, col=i_col_coords,
        par=i, mult=torch_tensor(i_mult[j], device=device))))
    }
  }

  torch_matrices <- list(
    A=torch_tensor(model$num_matrices[["A"]], device=device, dtype=dtypes[['A']]),
    Fm= torch_tensor(model$num_matrices[["Fm"]], device=device, dtype=dtypes[['Fm']]),
    S=torch_tensor(model$num_matrices[["S"]], device=device, dtype=dtypes[['S']]),
    Sk=torch_tensor(model$num_matrices[["Sk"]], device=device, dtype=dtypes[['Sk']]),
    K=torch_tensor(model$num_matrices[["K1_ref"]]+1-1, device=device, dtype=dtypes[['K']]),
    diag_n_p=torch_tensor(torch_diagflat(rep(1, model$meta_data$n_phenotypes + model$meta_data$n_latent)), device=device, dtype=dtypes[['diag_n_p']])
  )
  param_list <- list(A=NULL, Fm=NULL, S=NULL, Sk=NULL, K=NULL)
  .par_list <- list(A=NULL, Fm=NULL, S=NULL, Sk=NULL, K=NULL)

  torch_maps <- list(A=list(), Fm= list(), S=list(), Sk=list(), K=list())
  torch_masks <- list(
    A=torch_ones_like(torch_matrices[["A"]], dtype=dtypes[['A']], device=device),
    Fm= torch_ones_like(torch_matrices[["Fm"]], dtype=dtypes[['Fm']], device=device),
    S=torch_ones_like(torch_matrices[["S"]], dtype=dtypes[['S']], device=device),
    Sk=torch_ones_like(torch_matrices[["Sk"]], dtype=dtypes[['Sk']], device=device),
    K=torch_ones_like(torch_matrices[["K"]], dtype=dtypes[['K']], device=device)
  )
  for (i in torch_coords) {
    new_mat <- torch_zeros_like(torch_matrices[[i$mat_name]], device=device, dtype=dtypes[[i$mat_name]])
    new_mat[i$row, i$col] <- 1 * i$mult
    torch_masks[[i$mat_name]][i$row, i$col] <- 0.0
    param_name <- gsub("-", "", model$named_matrices[[i$mat_name]][i$row, i$col])
    if (param_name %in% names(torch_maps[[i$mat_name]])) {
      torch_maps[[i$mat_name]][[param_name]] <- torch_maps[[i$mat_name]][[param_name]]  + new_mat
    } else {
      torch_maps[[i$mat_name]][[param_name]] <- new_mat
      param_list[[i$mat_name]] <- c(param_list[[i$mat_name]], param_name)
      .par_list[[i$mat_name]] <- c(.par_list[[i$mat_name]], model$param_values[model$param_names == param_name])
    }
  }
  n_p <- model$meta_data$n_phenotypes + model$meta_data$n_latent

  # 3D tensors defining locations of paramters that require grad
  for (i in c("A", "Fm", "S", 'Sk', 'K')) {
    torch_maps[[i]] <- if (length(torch_maps[[i]]) > 0) {torch_dstack(torch_maps[[i]])} else {torch_zeros_like(torch_matrices[[i]], device=device, dtype=dtypes[[i]])}
  }
  # Reshape to 3D is necessary as we sum over 3rd axis later on
  for (i in c("A", "Fm", "S", 'Sk', 'K')) {
    if (length(torch_maps[[i]]$shape) == 2) {
      shape <- as.numeric(torch_maps[[i]]$shape)
      torch_maps[[i]] <- torch_reshape(torch_maps[[i]], c(shape[1], shape[2], 1))
    }
  }
  K2 <- torch_ones_like(torch_matrices[['K']], device=device, dtype=dtypes[['K']])
  for (i in seq_len(n_p)) {
    # Copy kurtosis of latent factors to K2 matrix
    coords <- .nd_to_2d_idx(n_p, i, i, i, i)
    if (model$named_matrices$K[coords$x, coords$y] != "0") {
      # K2[coords$x, coords$y]  <- model$num_matrices$K[coords$x, coords$y]
      K2[coords$x, coords$y]  <- 3.0
    }
  }
  base_matrices <- list(
    A=torch_mul(torch_tensor(model$num_matrices[["A"]], device=device, dtype=dtypes[['A']]), torch_masks[['A']]),
    Fm=torch_mul(torch_tensor(model$num_matrices[["Fm"]], device=device, dtype=dtypes[['Fm']]), torch_masks[['Fm']]),
    S=torch_mul(torch_tensor(model$num_matrices[["S"]], device=device, dtype=dtypes[['S']]), torch_masks[['S']]),
    Sk=torch_mul(torch_tensor(model$num_matrices[["Sk"]], device=device, dtype=dtypes[['Sk']]), torch_masks[['Sk']]),
    K=torch_tensor(model$num_matrices[["K1_ref"]]+1-1, device=device, dtype=dtypes[['K']]),
    K2=K2,
    diag_n_p=torch_tensor(torch_diagflat(rep(1, model$meta_data$n_phenotypes + model$meta_data$n_latent)), device=device, dtype=dtypes[['diag_n_p']])
  )
  # Zeros can produce NAN gradients, therefore set values in S matrices to very low values
  base_matrices[['S']] <- (base_matrices[['S']] + torch_tensor(1e-16, device=device))
  # If the matrix (A, Fm, S, Sk, or K) does not have free parameters, use a constant 1 to ensure code will always work
  .par_list <- list(
    A=if (length(param_list[['A']]) > 0) {torch_tensor(as.numeric(.par_list[['A']]), requires_grad = TRUE, device=device, dtype=dtypes[['A']])} else {torch_tensor(1, device=device, dtype=dtypes[['A']])},
    Fm=if (length(param_list[['Fm']]) > 0) {torch_tensor(as.numeric(.par_list[['Fm']]), requires_grad = TRUE, device=device, dtype=dtypes[['Fm']])} else {torch_tensor(1, device=device, dtype=dtypes[['Fm']])},
    S=if (length(param_list[['S']]) > 0) {torch_tensor(as.numeric(.par_list[['S']]), requires_grad = TRUE, device=device, dtype=dtypes[['S']])} else {torch_tensor(1, device=device, dtype=dtypes[['S']])},
    Sk=if (length(param_list[['Sk']]) > 0) {torch_tensor(as.numeric(.par_list[['Sk']]), requires_grad = TRUE, device=device, dtype=dtypes[['Sk']])} else {torch_tensor(1, device=device, dtype=dtypes[['Sk']])},
    K=if (length(param_list[['K']]) > 0) {torch_tensor(as.numeric(.par_list[['K']]), requires_grad = TRUE, device=device, dtype=dtypes[['K']])} else {torch_tensor(1, device=device, dtype=dtypes[['K']])}
  )
  torch_bounds <- list(L=list(), U=list())
  for (i in names(.par_list)) {
    if (.par_list[[i]]$requires_grad) {
      torch_bounds[["L"]][[i]] <- torch_tensor(as.numeric(model$bounds["L", param_list[[i]]]), device=device, dtype=dtypes[['bounds']])
      torch_bounds[["U"]][[i]] <- torch_tensor(as.numeric(model$bounds["U", param_list[[i]]]), device=device, dtype=dtypes[['bounds']])
    } else {
      # Non-grad parameters also need bounds in order for the shapes to match
      torch_bounds[["L"]][[i]] <- (torch_ones_like(.par_list[[i]], device=device, dtype=dtypes[['bounds']]) - 100)
      torch_bounds[["U"]][[i]] <- (torch_tensor(.par_list[[i]], device=device, dtype=dtypes[['bounds']]) + 100)
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
