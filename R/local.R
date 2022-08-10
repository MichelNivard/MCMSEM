MCMSEMversion <- "0.7.2"

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

# Turn ND (>2D) index into 2D index for use with M3 and M4
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