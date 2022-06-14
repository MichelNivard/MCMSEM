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
t4crossprod <- function(x){

  c(.m2m2v(tcrossprod(x)),.m3m2v(t(x%o%x%x%x)),.m4m2v(t(x%o%x%x%x%x%x)))

}

######## Compute Jacobian:
jac.fn <- function(par,model){

  n_p <- model$meta_data$n_phenotypes + model$meta_data$n_confounding
  # Model function
  ### Assign new parameter values to the matrices
  model$param_values <- par
  for (i in 1:length(model$param_coords)) {
    model$num_matrices[[model$param_coords[[i]][[1]]]][model$param_coords[[i]][[2]]] <- model$param_values[i] * model$param_coords[[i]][[3]]
  }
  # Extract matrices
  A  <- model$num_matrices[["A"]]
  Fm <- model$num_matrices[["Fm"]]
  S  <- model$num_matrices[["S"]]
  Sk <- model$num_matrices[["Sk"]]
  K <- model$num_matrices[["K"]]

  K[,] <- 0
  # there are some non 0 entries in S4, fix those using existing K1_ref
  K[model$num_matrices[["K1_ref"]]] <- 1
  # these are function of S2 matrix
  K <- sqrt(S) %*% K %*% (sqrt(S) %x% sqrt(S) %x% sqrt(S))
  for (i in 1:model$meta_data$n_confounding) {
    K[i, i + (i-1)*(n_p) + (i-1)*((n_p)^2)]  <- 3
  }
  # Re-enter values for K
  for (i in 1:length(model$param_coords)) {
    if (model$param_coords[[i]][[1]] == "K") {
      K[model$param_coords[[i]][[2]]] <- model$param_values[i]
    }
  }

  ###### Compute the observed cov, cosk, and cokurt matrices #################
  ############################################### (see section 2.2 paper) ####
  M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
  M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))

  ### Loss function
  value <- c(.m2m2v(M2),.m3m2v(M3),.m4m2v(M4))
  value
}

### pull it together to make std errors:
std.err <- function(data,par,model){

  n <- nrow(data)
  # if there is too much data for spoeedly opperation, sample 250000 observations to base this on
  if(n > 250000){

    samp <- sample(1:n,250000,F)
    data <- data[samp,]
  }

  # observed cov between pseudo obsertvations ovver n-1 gets us cov betwene moments moments
  S.m <- cov(t(apply(scale(data,center = T,scale = F),1,t4crossprod)))/(n-1)
  # weights matrix, works well if N is large
  W <- solve(S.m)

  G <- jacobian(func = jac.fn,x = par,model=model)
  Asycov <- solve(t(G)%*%W%*%G)

  se <- sqrt(diag(Asycov))

  se

}



