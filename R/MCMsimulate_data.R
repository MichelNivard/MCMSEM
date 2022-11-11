##################################################
# DATA SIMULATION
##################################################

simulate_data <- function(n=500000,
                          a=matrix(c(0.35, 0.35), ncol=1), b=matrix(c(1,-.1, .3,1), 2, 2),
                          shape=c(4, 0), df=c(0, 10),
                          asdataframe=FALSE, ...) {
  # a is an n_variable * n_latent matrix containing a parameters
  # b is an n_variable * n_variable matrix containing b parameters with 1 on diagonal:
  #  1   b2_1  Where b1_2 is the effect of x1 on x2, and b2_1 the effect of x2 on x1
  # b1_2  1
  # shape is a vector of shapes to be used for skewness (passed to rgamma), 0 is interpreted as no skewness
  # df is a vector of dfs to be used for kurtosis (passed to rt), 0 is interpreted as no kurtosis
  cat("  Note this data simulation is not exact and may contain sampling error.\n")
  if (n < 100000) {
    warning("Sampling error is likely worsened by the small sample size")
  }
  n_latent <- ncol(a)
  n_var <- nrow(a)
  if (ncol(b) != nrow(b))
    stop("b should be a square matrix")
  if (ncol(b) != n_var)
    stop("number of rows in b matrix should be qual to number of rows of a matrix")
  if (length(shape) != n_var)
    stop("shape should have a length equal to number of variables")
  if (length(df) != n_var)
    stop("df should have a length equal to number of variables")

  # Normal confounders for each x
  latent_f <- list()
  for (i in seq_len(n_latent)) {
    latent_f[[i]] <- rnorm(n)
  }
  latent_f <- do.call(cbind, latent_f)

  # Make variables:
  x <- list()
  for (i in seq_len(n_var)) {
     x[[i]] <- rowSums(a[i, ] * latent_f)
    if (shape[i] != 0)
      x[[i]] <- x[[i]] + scale(rgamma(n,shape=shape[i],scale=2)) # Non normal errors
    if (df[i])
      x[[i]] <- x[[i]] + scale(rt(n,df=df[i])) # Non normal errors
  }

  x <- do.call(cbind, x)

  # Let them concurrently influence each other
  if (asdataframe) {
    return(t(b %*% t(x)))
  } else {
    return(MCMdatasummary(t(b %*% t(x)), ...))
  }
}
