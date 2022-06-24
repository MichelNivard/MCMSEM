##################################################
# DATA SIMULATION
##################################################

simulate_data <- function(n=500000, a1=0.35, b1=0.3, b2=-.1, shape=4, df=10) {
  cat("  Note this data simulation is not exact and may contain sampling error.\n")
  if (n < 100000) {
    warning("Sampling error is likely worsened by the small sample size")
  }
  # Normal confounder
  f <- rnorm(n)

  # Non normal errors
  e1 <- scale(rgamma(n,shape=shape,scale=2))
  e2 <- scale(rt(n,df=df))

  # Make variables X1 and X2:
  x1 <- a1*f + e1
  x2 <- a1*f + e2

  # Let them concurrently influence each other
  hold <-matrix(c(1,b2, b1,1), 2, 2,byrow=T) %*% t(cbind(x1,x2))
  x1 <- hold[1,]
  x2 <- hold[2,]
  data <- cbind(x1,x2)
  return(data)
}