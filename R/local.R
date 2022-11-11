MCMSEMversion <- "0.23.0"

# Implemented loss functions
.get_lossfunc <- function(loss_type) {
  # For future reference, should we/someone ever want to create a custom loss function...
  # Note that this function returns a function!, which takes as input 2 tensors and calculates a 0D tensor contining the loss from those...
  #  in short:
  # x <- .get_lossfunc("customloss")
  # loss <- x(t1, t2)
  lossfuncs <- list(
    mse=nn_mse_loss(reduction='sum'),
    smooth_l1=nn_smooth_l1_loss(reduction='sum'),
    l1=nn_l1_loss(reduction='sum')
  )
  if (!(loss_type %in% names(lossfuncs))) {
    stop("loss_type should be one of c('mse', 'smooth_l1')")
  }
  return(lossfuncs[[loss_type]])
}

# Implemented optimizers for optimizer 1
.get_optimfunc <- function(optimizer) {
  # For future reference, should we/someone ever want to create a custom optimizer...
  # Similar to .get_lossfunc, this function returns a function, though a more complicated one, see https://torch.mlverse.org/docs/reference/optimizer.html
  # From that example, at the time of writing, the object optim_sgd2 would have to be returned from .get_optimfunc, i.e.
  # x <- .get_optimfunc("sgd2")
  # opt <- x(params, learning_rate=0.1)
  optimfuncs <- list(
    rprop=optim_rprop,
    sgd=optim_sgd,
    rmsprop=optim_rmsprop,
    asgd=optim_asgd,
    adam=optim_adam,
    adagrad=optim_adagrad,
    adadelta=optim_adadelta
  )
  if (!(optimizer %in% names(optimfuncs))) {
    stop("optimizer should be one of c(",paste0("'", names(optimfuncs), "'", collapse=", "),")")
  }
  return(optimfuncs[[optimizer]])
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

# Turn ND (>=2D) index into 2D index for use with M3 and M4
.nd_to_2d_idx <- function(nrows, x, y, ...) {
  # nrow must be provided as that determines stepsizes
  # x: first dim idx (i.e. row); y = second dim idx (i.e. col), ... = additional dim(s) idx(s)
  dim_idx <- list(...)
  for (i in seq_along(dim_idx)) {
    y <- y + (nrows^(i))*(dim_idx[[i]]-1)
  }
  return(list(x=x, y=y))
}

# The inverse of .nd_to_2d_idx
.twod_to_nd_idx <- function(nrows, x, y, ndims) {
  nd_coords <- list()
  nd_coords[[1]] <- x
  for (i in ndims:2) {
    nd_coords[[i]] <- floor((y - 1) / nrows^(i-2)) + 1
    y <- ((y - 1) %% nrows^(i-2)) + 1
  }
  return(nd_coords)
}

# Generate M2, M3, M4 comoment matrices
.get_comoments <- function(data, weights=NULL, debug=FALSE) {
  if (is.null(weights)) {
    if (debug) {cat(" - M2\n")}
    M2 <- cov(data)
    if (debug) {cat(" - M3\n")}
    M3 <- M3.MM(data)
    if (debug) {cat(" - M4\n")}
    M4 <- M4.MM(data)
  } else {
    # Custom implementation of weighted M2, M3, and M4 matrices
    weights <- torch_reshape(torch_tensor(weights), c(nrow(data),1))
    tens <- torch_tensor(data)

    weighted_mean <- torch_sum((tens*weights)/torch_sum(weights), dim=1)
    centred_tens <- weights*(tens - torch_reshape(weighted_mean, c(1, ncol(data))))
    if (debug) {cat(" - wM2\n")}
    M2 <- torch_matmul((centred_tens)$t(), centred_tens)/(torch_sum(weights))
    if (debug) {cat(" - wM3\n")}
    M3 <- torch_zeros(c(ncol(data), ncol(data)^2))
    for (i in seq_len(ncol(data))) {
      centred_tens_ <- centred_tens * torch_reshape(centred_tens[, i], c(nrow(data), 1))
      M3[, ((i-1)*ncol(data)+1):(i*ncol(data))] <- torch_matmul(centred_tens_$t(), centred_tens)
    }
    M3 <- M3/torch_sum(weights)
    if (debug) {cat(" - wM4\n")}
    M4 <- torch_zeros(c(ncol(data), ncol(data)^3))
    iter <- 0
    for (i in seq_len(ncol(data))) {
      for (j in seq_len(ncol(data))) {
        iter <- iter + 1
        centred_tens_ <- centred_tens * torch_reshape(centred_tens[, i], c(nrow(data), 1)) * torch_reshape(centred_tens[, j], c(nrow(data), 1))
        M4[, ((iter-1)*ncol(data)+1):(iter*ncol(data))] <- torch_matmul(centred_tens_$t(), centred_tens)
      }
    }
    M4 <- M4/torch_sum(weights)
  }
  return(list(M2=M2, M3=M3, M4=M4))
}

.twod_to_nd_idx <- function(nrows, x, y, ndims) {
  nd_coords <- list()
  nd_coords[[1]] <- x
  for (i in ndims:2) {
    nd_coords[[i]] <- floor((y - 1) / nrows^(i-2)) + 1
    y <- ((y - 1) %% nrows^(i-2)) + 1
  }
  return(nd_coords)
}

.torch_slowneckerproduct <- function(x, y, kronrow) {
  # This function is superceeded by .jit_slownecker$fn, but I'm leaving this here so you can read what happens in R language (note it is conensed, so copy-paste somewhere else and use newlines to make it more readable) should anything ever need changing and the person changing it is not comfortable in PyTorch
  device <- x$device
  out <- torch_tensor(matrix(0, nrow=nrow(x), ncol=sum(kronrow)), device=device)
  firstprod <- .torch_kron(y, y) # In case of memory issues: drop this line
  ncol <- 1
  for (j in seq_len(ncol(y)^3)) {
    if (kronrow[j]) {
      idx1 <- (j-1) %% ncol(y) + 1
      idx2 <- (floor((j-1)/ ncol(y)) %% ncol(y)) + 1
      idx3 <- (floor((j-1)/ ncol(y)^2) %% ncol(y)) + 1
      kroncol <- .torch_kron(firstprod[,idx1+idx2*ncol(y)], y[,idx3]) # In case of memory issues replace with: kroncol <- .torch_kron(.torch_kron(y[, idx1], y[, idx2]), y[, idx3])
      for (i in seq_len(nrow(out))) {
        out[i, ncol] <- torch_matmul(kroncol, x[i, ])
      }
      ncol <- ncol + 1
    }
  }
  return(out, device=device)
}
