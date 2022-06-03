# Fit MCM model
# Objective function
.objective <- function(.par, model, M2.obs, M3.obs, M4.obs){
  n_p <- model$meta_data$n_phenotypes + model$meta_data$n_confounding
  # Model function
  ### Assign new parameter values to the matrices
  model$param_values <- .par
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
  value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  return(value)
}


MCMfit <- function(model, data, compute_se=TRUE, bootstrap_type='two-step', bootstrap_iter=200,bootstrap_chunks=1000) {
  #TODO: Add arguments for fitting either x->y or y->x path as opposed to both (which should remain the default)
  #TODO: Expand manual
  #if (!(confounding %in% c('positive', 'negative', 'both')))
  #  stop("confounding should be one of c('positive', 'negative', 'both')")
  if (!(bootstrap_type %in% c('two-step', 'one-step')))
    stop("confounding should be one of c('two-step', 'one-step')")
  if (ncol(data) != 2)
    warning("Use of a dataframe with more than 2 columns is still experimental.")
  if (nrow(data) < 1000)
    stop("Currently only a dataframe with at least 1000 rows is supported.")
  if (ncol(data) != model$meta_data$n_phenotypes) {
    msg <- paste0("Model expected ", model$meta_data$n_phenotypes, " phenotypes but ", ncol(data), " were found in the data. Please create a new model with this dataset.")
    stop(msg)
  }
  #TODO: There is no way to do both currently
  #if (confounding == 'both') {
  #  model2 <- model$copy()  # Make a copy of the model instance first, as parameter values are modified inplace
  #  result_positive <- MCMSEM(model, data, confounding="positive", compute_se=compute_se, bootstrap_type=bootstrap_type, bootstrap_iter=bootstrap_iter,bootstrap_chunks=bootstrap_chunks)
  #  result_negative <- MCMSEM(model2, data, confounding="negative", compute_se=compute_se, bootstrap_type=bootstrap_type, bootstrap_iter=bootstrap_iter,bootstrap_chunks=bootstrap_chunks)
  #  return(list(positive_confounder=result_positive, negative_confounder=result_negative))
  #}
  # Scale data
  if (model$meta_data$scale_data) {
      for (i in 1:ncol(data)) {
        data[,i] <- scale(data[,i])
      }
  }

  if (compute_se)
    model_copy <- model$copy()

  # Obtain covariance, coskewness and cokurtosis matrices
  M2.obs <- cov(data)
  M3.obs <- M3.MM(data)
  M4.obs <- M4.MM(data)

  # Specify starting values
  start <- model$param_values

  # Specify upper and lower bound of parameters
  L <- as.numeric(model$bounds["L", ])
  U <- as.numeric(model$bounds["U", ])

  # Obtain estimates with optimizer
  nlminb.out <-nlminb(start,objective = .objective,model=model,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,lower = L, upper = U)

  # Store estimates including minimization objective, using this to evaluate/compare fit
  results        <-  as.data.frame(matrix(c(nlminb.out$par, nlminb.out$objective), nrow = 1))

  ## NOTE:
  # When fitting to real data, we compare model with pos. confounder with
  # model with neg. confounder (see matrix A; section 2.2 paper). That is,
  # the sign in front of one of the a1s in A becomes (-) (see code above).
  # We run both models and choose the model with the lowest minimize.obj
  # This might be something that should be built in automatically.

  if (compute_se) {
    # Matrix where bootstraps will be stored
    pars.boot <- matrix(NA,bootstrap_iter,length(model$param_values))
    # Lower and Upper bounds
    L <- as.numeric(model$bounds["L", ])
    U <- as.numeric(model$bounds["U", ])
    cat("Starting bootstrap MCMSEM\n")
    pb <- txtProgressBar(0, bootstrap_iter, style = 3, width=min(c(options()$width, 107)))
    if (bootstrap_type == 'one-step') {
      #### BOOT 1:  NORMAL BOOTSTRAP
      ##############################
      timestart <- Sys.time()

      # Matrix where bootstraps will be stored
      # Bootstrap
      for (i in 1:bootstrap_iter){
        ## Progress bar stuff
        setTxtProgressBar(pb, i)
        model2 <- model_copy$copy()  # Create new empty model
        #1. Sample from data with replacement
        boot <-   sample(1:nrow(data),nrow(data),T)
        sample <- data[boot,]

        #2. Get covariance, coskewness and cokurtosis matrices
        M2.obs <-   cov(sample)
        M3.obs <- M3.MM(sample)
        M4.obs <- M4.MM(sample)

        #3. Fit model
        # Estimate parameters with model function specified above
        nlminb.out <-nlminb(model2$param_values,objective = .objective,model=model2,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,lower = L, upper = U)
        # Store point estimates of bootstraps
        pars.boot[i,] <- nlminb.out$par

      }
      close(pb)
      SEs_boot <- apply(pars.boot, 2, sd)
      timeend <- Sys.time()
      boot1 <- timeend-timestart
    } else if (bootstrap_type == 'two-step') {
      #### BOOT 2: TWO STEP BOOTSTRAP
      ##############################
      starttime <- Sys.time()

      ### STEP 1
      # 1. Bind the data to a random sample binning people into "bootstrap_chunks" groups
      step1 <- cbind(data,sample(1:bootstrap_chunks,nrow(data),replace=T))
      colnames(step1)[ncol(step1)] <- "group"
      step1 <- as.data.frame(step1)

      # 2. Get covariance, coskenwess and cokurtosis matrices per group
      sample.cov  <- aggregate(1:nrow(step1), by=list(step1$group), function(s) cov(step1[s, colnames(data)]))
      sample.cosk <- aggregate(1:nrow(step1), by=list(step1$group), function(s) M3.MM(as.matrix(step1[s, colnames(data)])))
      sample.cokr <- aggregate(1:nrow(step1), by=list(step1$group), function(s) M4.MM(as.matrix(step1[s, colnames(data)])))
      pars.boot2 <- matrix(NA,nrow=bootstrap_iter,ncol=length(model$param_values))

      ### STEP 2

      for (i in 1:bootstrap_iter){
        ## Progress bar stuff
        setTxtProgressBar(pb, i)
        model2 <- model_copy$copy()  # Create new empty model
        # 3. Sample cov/cosk/cokrt matrices and obtain mean of the sampled matrices
        #to use as cov/cosk/cokrt matrix in the model
        boot <-   sample(1:bootstrap_chunks,bootstrap_chunks,T)
        M2.obs <-   matrix(colMeans(sample.cov [boot,-1]),ncol(data),ncol(data), byrow=T)
        M3.obs <-   matrix(colMeans(sample.cosk[boot,-1]),ncol(data),ncol(data)^2, byrow=T)
        M4.obs <-   matrix(colMeans(sample.cokr[boot,-1]),ncol(data),ncol(data)^3, byrow=T)


        # 4. Fit model,
        nlminb.out <-nlminb(model2$param_values,objective = .objective,model=model,M2.obs=M2.obs,M3.obs=M3.obs,M4.obs=M4.obs,lower = L, upper = U)
        pars.boot2[i,] <- nlminb.out$par

      }
      close(pb)
      SEs_boot <- apply(pars.boot2, 2, sd)
      endtime <- Sys.time()
      boot2 <- endtime - starttime
    }
    #cat("\n")  # Prevent things from printing over completed progress bar
    # Table of point estimates and SE's
    results <- rbind(results, c(SEs_boot,NA))
  }

  colnames(results) <- c(model$param_names, "mimize.obj")
  rownames(results) <- if(compute_se) c("est", "se") else c("est")
  return(results)
}

