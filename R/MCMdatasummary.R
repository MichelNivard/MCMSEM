MCMdatasummary <- function(data=NULL, path=NULL, scale_data=TRUE, prep_asymptotic_se=TRUE, weights=NULL) {
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
      data_was_scaled <- all(round(apply(data, 2, mean), 1) == 0) & all(round(apply(data, 2, sd), 1) == 1)
      if (!(data_was_scaled)) {
        data <- apply(data, 2, scale)
      }
    } else {
      data_was_scaled <- TRUE
    }
    data <- as.matrix(data)
    weighted <- !(is.null(weights))
    if (weighted) {weightsum <- sum(weights)} else {weightsum <- nrow(data)}
    comoments <- get_comoments(data, weight=weights)
    n <- nrow(data)

    if (prep_asymptotic_se) {
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
      # S.m <- cov(t(S.m))/(n-1)
      # Replace this with torch_cov(S.m) / (n - 1)  once that is implemented in torch for R
      S.m <- apply(scale(data, center = T, scale = F),1, .t4crossprod,
                   idx=c(dim2locs, length(as.vector(dim2)) + dim3locs,  length(as.vector(dim2)) + length(as.vector(dim3)) + dim4locs),
                   use_skewness=TRUE, use_kurtosis=TRUE)
      S.m <- torch_transpose(torch_tensor(S.m, device=torch_device("cpu")), 1, 2)
      S.m <- torch_subtract(S.m, torch_reshape(torch_mean(S.m, dim=1), c(1, S.m$shape[2])))
      S.m <- torch_matmul(torch_transpose(S.m, 1, 2), S.m)  / (S.m$shape[1] - 1) / (n - 1)
      idx <- list(
        idx=1:ncol(S.m),
        idx_nokurt=1:(ncol(S.m)-length(dim4locs)),
        idx_noskew=c(1:length(dim2locs), (length(dim2locs)+length(dim3locs)+1):ncol(S.m)),
        idx_nokurt_noskew=1:length(dim2locs)
      )
      SE <- list(
        computed=TRUE,
        S.m=as.matrix(S.m),
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




