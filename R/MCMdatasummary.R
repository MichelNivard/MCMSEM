MCMdatasummary <- function(data=NULL, path=NULL, weights=NULL, scale_data=TRUE, prep_asymptotic_se=TRUE, use_skewness=TRUE,
                           use_kurtosis=TRUE, debug=FALSE, low_memory=FALSE) {
  if (is.null(data) & is.null(path)) {
    stop("Either argument data (to generate a new summary) or argument path (to load an existing summary) should be provided")
  }
  if (is.null(data)) {
    res <- mcmdataclass()
    res$load(path)
    return(res)
  } else if (is.null(path) & is.character(data) & (length(data) == 1)) {
    warning("A string of length 1 was provided to the argument data, I'm assuming this is a path to an existing summary file.")
    res <- mcmdataclass()
    res$load(data)
    return(res)
  } else {
    if (nrow(data) < 1000)
      stop("Currently only a dataframe with at least 1000 rows is supported.")
    if (ncol(data) < 2)
      stop("At least two columns in data are required.")
    if (any(apply(data, 2, is.character)))
      stop("Numeric column(s) found:",paste0(names(which(apply(data, 2, is.character))), collapse=", "), "\n   Data must only contain numeric columns")
    if (any(round(rowSums(abs(cor(data))), 1) == 1.0)) {
      idcol <- which(round(rowSums(abs(cor(data))), 1) == 1.0)
      stop(paste0("It seems like ",names(idcol), " is an ID column. \n  Consider using ", deparse(substitute(data)), '[,-', unname(idcol),']'))
    }
    if (!is.null(weights)) {if (length(weights) != nrow(data)) {stop("Number of weights should be equal to the number of rows in the dataframe")}}
    ncol <- ncol(data)
    colnames <- if (is.data.frame(data)) {colnames(data)} else {paste0("x", seq_len(ncol(data)))}
    if (is.matrix(data) & !(is.null(colnames(data)))) {colnames <- colnames(data)}

    # Scale data
    if (scale_data) {
      if (debug) {cat("MCMdatasummary Scaling data\n")}
      data_was_scaled <- all(round(apply(data, 2, mean), 1) == 0) & all(round(apply(data, 2, sd), 1) == 1)
      if (!(data_was_scaled)) {
        data <- apply(data, 2, scale)
      }
    } else {
      data_was_scaled <- TRUE
    }
    if (debug) {cat("MCMdatasummary Converting to matrix\n")}
    data <- as.matrix(data)
    weighted <- !(is.null(weights))
    if (weighted) {weightsum <- sum(weights)} else {weightsum <- nrow(data)}
    if (debug) {cat("MCMdatasummary calculating comoments:\n")}
    comoments <- .get_comoments(data, weight=weights, debug=debug)
    n <- nrow(data)

    if (prep_asymptotic_se) {
      if (debug) {cat("MCMdatasummary calculating S.m:\n")}
      # if there is too much data for spoeedly opperation, sample 100000 observations to base this on
      if(n > 100000){
        if (debug) {cat(" - resampling\n")}
        samp <- sample(1:n,100000,F)
        data <- torch_tensor(data[samp, ])
      } else {
        data <- torch_tensor(data)
      }
      if (debug) {cat(" - calculating dimlocs\n")}
      # observed cov between pseudo obsertvations ovver n-1 gets us cov betwene moments moments
      dim2locs <- torch_tensor(.dimlocations(ncol(data), dims=2) - 1, dtype=torch_int64())
      dim3locs <- torch_tensor(.dimlocations(ncol(data), dims=3) - 1, dtype=torch_int64())
      dim4locs <- torch_tensor(.dimlocations(ncol(data), dims=4) - 1, dtype=torch_int64())

      # R style: t4crossprod
      #x2 <- x %o% x %x% x %x% x
      #x3 <- x2 %x% x
      #return(c(as.vector(t(x2)), as.vector(t(x3)), as.vector(t(x4)))[idx])
      # apply that to scaled data across rows
      # with idx = c(as.numeric(dim2locs), length(as.vector(dim2)) + as.numeric(dim3locs),  length(as.vector(dim2)) + length(as.vector(dim3)) + as.numeric(dim4locs)) + 1
      if (debug) {cat(" - calculating t4crossprod\n")}
      .jit_t4crossprod <- jit_compile(.jit_funcs[['t4crossprod']])
      if (use_kurtosis & use_skewness) {
        S.m <- .jit_t4crossprod$fn(data, dim2locs, dim3locs, dim4locs)
      } else if (use_kurtosis) {
        S.m <- .jit_t4crossprod$fn(data, dim2locs, torch_tensor(0), dim4locs)
      } else if (use_skewness) {
        S.m <- .jit_t4crossprod$fn(data, dim2locs, dim3locs, torch_tensor(0))
      } else {
        S.m <- .jit_t4crossprod$fn(data, dim2locs, torch_tensor(0), torch_tensor(0))
      }

      # S.m <- cov(t(S.m))/(n-1)
      # Replace this with torch_cov(S.m) / (n - 1)  once that is implemented in torch for R
      if (low_memory) {
        if (debug) {cat(" - calculating S.m covariance using lowmem\n")}
        # This currently takes forever
        #TODO: chunked covariance calculation
        covlowmem <- jit_compile(.jit_funcs[['covlowmem']])$fn
        S.m <- covlowmem(S.m) / (n - 1)
      } else {
        if (debug) {cat(" - calculating S.m covariance\n")}
        S.m <- torch_subtract(S.m, torch_reshape(torch_mean(S.m, dim=1), c(1, S.m$shape[2])))
        S.m <- torch_matmul(torch_transpose(S.m, 1, 2), S.m)  / (S.m$shape[1] - 1) / (n - 1)
      }
      if (debug) {cat(" - assigning S.m indices\n")}
      idx <- list()
      if (use_skewness & use_kurtosis) {
        idx[['idx']] <- 1:ncol(S.m)
        idx[['idx_nokurt']] <- 1:(ncol(S.m)-length(dim4locs))
        idx[['idx_noskew']] <- c(1:length(dim2locs), (length(dim2locs)+length(dim3locs)+1):ncol(S.m))
        idx[['idx_nokurt_noskew']] <- 1:length(dim2locs)
      } else if (use_kurtosis) {
        idx[['idx_noskew']] <- 1:ncol(S.m)
        idx[['idx_nokurt_noskew']] <- 1:length(dim2locs)
      } else if (use_skewness) {
        idx[['idx_nokurt']] <- 1:ncol(S.m)
        idx[['idx_nokurt_noskew']] <- 1:length(dim2locs)
      } else {
        idx[['idx_nokurt_noskew']] <- 1:length(dim2locs)
      }
      if (debug) {cat(" - converting S.m to R matrix\n")}
      SE <- list(
        computed=TRUE,
        S.m=S.m,
        idx=idx
      )
    } else {
      SE <- list(computed=FALSE, S.m=NULL, idx=NULL)
    }

    return(mcmdataclass(meta_data=list(scale_data=scale_data, data_was_scaled=data_was_scaled, weighted=weighted, ncol=ncol, colnames=colnames, N=n, weightsum=weightsum),
                        M2=comoments$M2, M3=comoments$M3, M4=comoments$M4, SE=SE))
  }
}

MCMsavesummary <- function(summaryobj, path) {
  summaryobj$save(path)
}




