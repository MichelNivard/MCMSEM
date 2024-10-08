.gen_matrices <- function(par, n_p, n_f, base_value=0) {
  matrices <- list()
  ############## A
  A <- matrix(rep(base_value, (n_p+n_f)^2), n_p+n_f)
  for (i in seq_along(par[['a']])) {
    a_val <- par[['a']][i]
    A[((i-1) %% n_p)+n_f+1, floor((i-1)/n_p)+1] <- a_val
  }
  B <- matrix(rep(base_value, (n_p)^2), n_p)
  iter <- 1
  for (i in seq_len(n_p)) {
    for (j in seq_len(n_p)) {
      if (i != j) {
        B[i, j] <- par[['b']][iter]
        iter <- iter + 1
      }
    }
  }
  A[(n_f+1):nrow(A), (n_f+1):ncol(A)] <- B
  matrices[['A']] <- A
  ############## Fm
  Fm <- matrix(rep(base_value, n_p * (n_p+n_f)), n_p)
  for (i in seq_len(n_p)) {
    Fm[i, i+n_f] <- if (is.character(base_value)) "1" else 1
  }
  matrices[['Fm']] <- Fm
  ############## S
  S <- matrix(rep(base_value, (n_p+n_f)^2), n_p+n_f)
  for (i in seq_len(n_f)) {
    S[i, i] <- if (is.character(base_value)) "1" else 1
  }
  for (i in seq_len(n_p)) {
    S[i+n_f, i+n_f] <- par[['s']][i]
  }
  matrices[['S']] <- S
  ############## Sk
  Sk <- matrix(rep(base_value, (n_p+n_f)^3), n_p + n_f)
  for (i in seq_len(n_p)) {
    Sk[i+n_f, i+n_f+(i-1+n_f)*(n_p+n_f)]  <- par[['sk']][i]
  }
  matrices[['Sk']] <- Sk
  ############## K
  K <- matrix(rep(base_value, (n_p+n_f)^4), n_p + n_f)
  if (!(is.character(base_value))) {
    K1_ref <- rnorm(2000)
    for (i in seq_len(n_p + n_f - 1)) {
      K1_ref <- cbind(K1_ref, rnorm(2000))
    }
    # K <- sqrt(S) %*% K %*% (sqrt(S) %x% sqrt(S) %x% sqrt(S))
    # Note this following solution only works if S is diagonal, assuming that is the case here:
    skron <- diag(S) %x% diag(S) %x% diag(S)
    for (i in seq_len(nrow(K))) {
      K[i, ] <- K[i, ] * skron
    }

    matrices[['K1_ref']] <- round(M4.MM(K1_ref)) != 0
  }
  for (i in seq_len(n_f)) {
    coords <- .nd_to_2d_idx(n_p+n_f, i, i, i, i)
    K[coords$x, coords$y] <- if (is.character(base_value)) as.character(3 * as.numeric(S[i, i]) ^ 2) else 3 * S[i, i] ^ 2
  }
  for (i in seq_len(n_p)) {
    coords <- .nd_to_2d_idx(n_p+n_f, i+n_f, i+n_f, i+n_f, i+n_f)
    K[coords$x, coords$y]  <- par[['k']][i]
  }

  matrices[['K']] <- K
  return(matrices)
}