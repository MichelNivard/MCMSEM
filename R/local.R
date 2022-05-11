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



###  .fn is only for reference in this branch.
.fn <- function(par,M2.obs,M3.obs,M4.obs,confounding){
  # Model function
  ### Assign the parameters
  a1  <- par[1]
  b1  <- par[2]
  b2  <- par[3]
  s1  <- par[4]
  s2  <- par[5]
  sk1 <- par[6]
  sk2 <- par[7]
  k1  <- par[8]
  k2  <- par[9]

  ####### Specify model matrices ##########################################
  ############### For explanation of each matrix, see section 2.2 of paper

  ###  A matrix (positive and negative conounder; turn-on only 1 !!)
  if (confounding == 'positive') {
  # positive
    A <- matrix(c(0,  0,  0,
                  a1, 0, b2,
                  a1, b1, 0), 3,3,byrow = T)
  } else if (confounding == 'negative') {
    A <- matrix(c(0,  0,  0,
                  -a1, 0, b2,
                  a1, b1, 0), 3,3,byrow = T)
  }

  ### F matrix
  Fm <- matrix(c(0, 1, 0,
                 0, 0, 1),  2,3,byrow = T)

  ###  S2 matrix (here S):
  S <- matrix(c(1, 0,  0,
                0, s1, 0,
                0, 0,  s2), 3,3,byrow = T)

  ###  S3 matrix (here Sk)
  Sk <- matrix(0,3,9)
  Sk[2,5] <- sk1
  Sk[3,9] <- sk2

  ###  S4 matrix (here K)
  K <- matrix(0,3,27)
  # there are some non 0 entries in S4
  K[ round(M4.MM(cbind(rnorm(2000), rnorm(2000), rnorm(2000)))) != 0] <- 1
  # these are function of S2 matrix
  K <- sqrt(S) %*% K %*% (sqrt(S) %x% sqrt(S) %x% sqrt(S))
  K[1,1] <- 3
  K[2,14] <- k1
  K[3,27] <- k2


  ###### Compute the observed cov, cosk, and cokurt matrices #################
  ############################################### (see section 2.2 paper) ####

  M2 <- Fm %*% solve(diag(3) - A) %*% S %*%  t(solve(diag(3)-A))  %*% t(Fm)
  M3 <- Fm %*% solve(diag(3) - A) %*% Sk %*% (t(solve(diag(3)-A)) %x% t(solve(diag(3)-A))) %*% (t(Fm) %x% t(Fm))
  M4 <- Fm %*% solve(diag(3) - A) %*% K %*%  (t(solve(diag(3)-A)) %x% t(solve(diag(3)-A))  %x%  t(solve(diag(3)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))

  ### Loss function
  value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)

}
