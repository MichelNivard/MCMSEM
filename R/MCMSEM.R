##################################################
# DATA SIMULATION
##################################################
simulate_data <- function(n=500000, a1=0.35, b1=0.3, b2=-.1, shape=4, df=10) {
  cat("  Note this data simulation is not exact and may contain sampling error.\n")
  if (n < 100000) {
    warning("Sampling error is likely worsened by the small sample size")
  }
  # Pick parameter values
  n     <- 500000
  a1    <- 0.35
  b1    <- 0.3
  b2    <- -0.1
  # shape <- 4   medium skew
  # df    <- 10  medium kurt

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


##################################################
# MODEL
##################################################
MCMSEM <- function(data, confounding="positive", compute_se=TRUE, bootstrap_type='two-step', bootstrap_iter=200) {
  if (!(confounding %in% c('positive', 'negative')))
    stop("confounding should be one of c('positive', 'negative')")
  if (!(bootstrap_type %in% c('two-step', 'one-step')))
    stop("confounding should be one of c('two-step', 'one-step')")
  if (ncol(data) != 2)
    stop("Currently only a dataframe with 2 columns is supported.")
  if (nrow(data) < 1000)
    stop("Currently only a dataframe with at least 1000 rows is supported.")
  # Obtain covariance, coskewness and cokurtosis matrices
  M2.obs <- cov(data)
  M3.obs <- M3.MM(data)
  M4.obs <- M4.MM(data)

  # Specify starting values
  start <- c(.2,.2,.2,1,1,M3.obs[1,1],M3.obs[2,4],M4.obs[1,1],M4.obs[2,8])

  # Specify upper and lower bound of parameters
  L <- c(-1,-.5,-.5,0.01,0.01,-5,-5,0,0)
  U <- c(1,1,1,2,2,18,18,100,100)

  # Obtain estimates with optimizer
  nlminb.out <-nlminb(start,objective = .fn,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,confounding=confounding,lower = L, upper = U)

  # Store estimates including minimization objective, using this to evaluate/compare fit
  results        <-  c(nlminb.out$par, nlminb.out$objective)

  ## NOTE:
  # When fitting to real data, we compare model with pos. confounder with
  # model with neg. confounder (see matrix A; section 2.2 paper). That is,
  # the sign in front of one of the a1s in A becomes (-) (see code above).
  # We run both models and choose the model with the lowest minimize.obj
  # This might be something that should be built in automatically.

  if (compute_se) {
    # Matrix where bootstraps will be stored
    pars.boot <- matrix(NA,bootstrap_iter,9)
    # Lower and Upper bounds
    L <- c(-1,-.5,-.5,0.01,0.01,-5,-5,0,0)
    U <- c(1,1,1,2,2,18,18,100,100)

    if (bootstrap_type == 'one-step') {
      #### BOOT 1:  NORMAL BOOTSTRAP
      ##############################
      timestart <- Sys.time()

      # Matrix where bootstraps will be stored

      # Bootstrap
      for (i in 1:bootstrap_iter){
        #1. Sample from data with replacement
        boot <-   sample(1:nrow(data),nrow(data),T)
        sample <- data[boot,]

        #2. Get covariance, coskewness and cokurtosis matrices
        M2.obs <-   cov(sample)
        M3.obs <- M3.MM(sample)
        M4.obs <- M4.MM(sample)

        #3. Fit model
        # Start values
        start <- c(.2,.2,.2,1,1,M3.obs[1,1],M3.obs[2,4],M4.obs[1,1],M4.obs[2,8])

        # Estimate parameters with model function specified above
        nlminb.out <-nlminb(start,objective = .fn,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,confounding=confounding,lower = L, upper = U)

        # Store point estimates of bootstraps
        pars.boot[i,] <- nlminb.out$par

      }
      SEs_boot <- apply(pars.boot, 2, sd)
      timeend <- Sys.time()
      boot1 <- timeend-timestart
    } else if (bootstrap_type == 'two-step') {
      #### BOOT 2: TWO STEP BOOTSTRAP
      ##############################
      starttime <- Sys.time()

      ### STEP 1
      # 1. Bind the data to a random sample binning people into 1000 groups
      step1 <- cbind(data,sample(1:1000,nrow(data),replace=T))
      colnames(step1)[3] <- "group"
      step1 <- as.data.frame(step1)

      # 2. Get covariance, coskenwess and cokurtosis matrices per group
      sample.cov  <- aggregate(1:nrow(step1), by=list(step1$group), function(s)   cov(matrix(c(step1$x1[s],step1$x2[s]),ncol=2,byrow=F)))
      sample.cosk <- aggregate(1:nrow(step1), by=list(step1$group), function(s) M3.MM(matrix(c(step1$x1[s],step1$x2[s]),ncol=2, byrow=F)))
      sample.cokr <- aggregate(1:nrow(step1), by=list(step1$group), function(s) M4.MM(matrix(c(step1$x1[s],step1$x2[s]),ncol=2, byrow=F)))
      pars.boot2 <- matrix(NA,bootstrap_iter,9)

      ### STEP 2

      for (i in 1:bootstrap_iter){

        # 3. Sample cov/cosk/cokrt matrices and obtain mean of the sampled matrices
        #to use as cov/cosk/cokrt matrix in the model
        boot <-   sample(1:1000,1000,T)
        M2.obs <-   matrix(colMeans(sample.cov [boot,-1]),2,2, byrow=T)
        M3.obs <-   matrix(colMeans(sample.cosk[boot,-1]),2,4, byrow=T)
        M4.obs <-   matrix(colMeans(sample.cokr[boot,-1]),2,8, byrow=T)


        # 4. Fit model,
        start <- c(.2,.2,.2,1,1,M3.obs[1,1],M3.obs[2,4],M4.obs[1,1],M4.obs[2,8])
        nlminb.out <-nlminb(start,objective = .fn,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,confounding=confounding,lower = L, upper = U)
        pars.boot2[i,] <- nlminb.out$par

      }
      SEs_boot <- apply(pars.boot2, 2, sd)
      endtime <- Sys.time()
      boot2 <- endtime - starttime
    }

    # Table of point estimates and SE's
    results <- rbind(results, c(SEs_boot,NA))
  }

  colnames(results) <- c("a", "b1", "b2", "vare1", "vare2", "skewe1", "skewe2", "kurte1", "kurte2", "mimize.obj")
  rownames(results) <- if(compute_se) c("est", "se") else c("est")
  return(results)
}

