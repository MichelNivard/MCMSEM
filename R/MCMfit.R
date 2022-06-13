# Fit MCM model
# Objective function
.torch_objective <- function(.par, model, torch_coords, M2.obs, M3.obs, M4.obs) {
  torch_A  <- torch_tensor(model$num_matrices[["A"]])
  torch_Fm <- torch_tensor(model$num_matrices[["Fm"]])
  torch_S  <- torch_tensor(model$num_matrices[["S"]])
  torch_Sk <- torch_tensor(model$num_matrices[["Sk"]])
  torch_K <- torch_tensor(model$num_matrices[["K"]])
  torch_matrices <- list(
    A=torch_A,
    Fm= torch_Fm,
    S=torch_S,
    Sk=torch_Sk,
    K=torch_K,
    diag_n_p=torch_diagflat(rep(1, model$meta_data$n_phenotypes + model$meta_data$n_confounding))
  )
  n_p <- model$meta_data$n_phenotypes + model$meta_data$n_confounding
  # Model function
  ### Assign new parameter values to the matrices
  for (i in torch_coords) {
    torch_matrices[[i$mat_name]][i$row, i$col] <- .par[i$par] * i$mult
  }
  # Extract matrices
  A  <- torch_matrices[["A"]]
  Fm <- torch_matrices[["Fm"]]
  S  <- torch_matrices[["S"]]
  Sk <- torch_matrices[["Sk"]]
  K <- torch_matrices[["K"]]
  diag_n_p <- torch_matrices[['diag_n_p']]

  K[,] <- 0
  # there are some non 0 entries in S4, fix those using existing K1_ref
  K[model$num_matrices[["K1_ref"]]] <- 1
  # these are function of S2 matrix
  sqrtS <- torch_sqrt(S)
  K <- torch_matmul(torch_matmul(sqrtS, K), .torch_kron(sqrtS, .torch_kron(sqrtS, sqrtS)))
  for (i in 1:model$meta_data$n_confounding) {
    K[i, i + (i-1)*(n_p) + (i-1)*((n_p)^2)]  <- 3
  }
  # Re-enter values for K
  for (i in torch_coords) {
    if (i$mat_name == "K") {
      K[i$row, i$col] <- .par[i$par] * i$mult
    }
  }

  ###### Compute the observed cov, cosk, and cokurt matrices #################
  ############################################### (see section 2.2 paper) ####
  diag_np_a_inv <- torch_inverse(diag_n_p - A)
  diag_np_a_inv_t <- torch_transpose(torch_inverse(diag_n_p - A), 1, 2)
  # M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), S), diag_np_a_inv_t), torch_transpose(Fm, 1, 2))
  # M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
  M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), Sk), .torch_kron(diag_np_a_inv_t, diag_np_a_inv_t)), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
  # M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))
  M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), K), torch_matmul(.torch_kron(.torch_kron(diag_np_a_inv_t, diag_np_a_inv_t), diag_np_a_inv_t), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))

  ### Loss function
  # value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  value <- torch_sum((torch_square(torch_tril(M2.obs) - torch_tril(M2))))  +  torch_sum((torch_square(torch_tril(M3.obs$reshape(c(2,2,2))) - torch_tril(M3$reshape(c(2,2,2)))))) + torch_sum((torch_square(torch_tril(M4.obs$reshape(c(2,2,2,2)))- torch_tril(M4$reshape(c(2,2,2,2))))))
  return(value)
}

.torch_fit <- function(model, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, max_iter, optim_tol) {
  .par <- torch_tensor(model$param_values, requires_grad=TRUE)
  optim <- optim_rprop(.par, lr=learning_rate)
  for (i in 1:max_iter) {
    optim$zero_grad()
    loss <- .torch_objective(.par, model, torch_coords, M2.obs, M3.obs, M4.obs)
    loss$backward()
    optim$step()
    #if (i > 1) cat(paste0("step",i,": loss", round(as.numeric(loss), 2), "diff:", round(as.numeric(loss-previous_loss), 2), "\n"))
    if (i > 2) {if (as.numeric(torch_less_equal(torch_abs(previous_loss-loss), optim_tol))) break}
    previous_loss <- loss
  }
  return(.par)
}

MCMfit <- function(model, data, compute_se=TRUE, se_type='asymptotic', bootstrap_iter=200,bootstrap_chunks=1000,
                   learning_rate=0.005, max_iter=5000, optim_tol=1e-8) {
  #TODO: Add arguments for fitting either x->y or y->x path as opposed to both (which should remain the default)
  #TODO: Expand manual
  #if (!(confounding %in% c('positive', 'negative', 'both')))
  #  stop("confounding should be one of c('positive', 'negative', 'both')")
  if (!(se_type %in% c('two-step', 'one-step','asymptotic')))
    stop("se_type should be one of c('two-step', 'one-step','asymptotic')")
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
  data <- as.matrix(data)
  # Scale data
  if (model$meta_data$scale_data) {
      for (i in 1:ncol(data)) {
        data[,i] <- scale(data[,i])
      }
  }

  if (compute_se) ###????
    model_copy <- model$copy()

  # Obtain covariance, coskewness and cokurtosis matrices
  M2.obs <- torch_tensor(cov(data))
  M3.obs <- torch_tensor(M3.MM(data))
  M4.obs <- torch_tensor(M4.MM(data))

  torch_coords <- .get_torch_coords(model)
  # Obtain estimates with optimizer
  .par <- .torch_fit(model, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, max_iter, optim_tol)
  # Store estimates including minimization objective, using this to evaluate/compare fit
  results <-  as.data.frame(matrix(as.numeric(.par), nrow = 1))

  ## NOTE:
  # When fitting to real data, we compare model with pos. confounder with
  # model with neg. confounder (see matrix A; section 2.2 paper). That is,
  # the sign in front of one of the a1s in A becomes (-) (see code above).
  # We run both models and choose the model with the lowest minimize.obj
  # This might be something that should be built in automatically.

  if (compute_se) {

    if(se_type == 'asymptotic'){
      SEs <- .std.err(data=data,par=as.numeric(.par),model=model)
    }

    if(se_type != 'asymptotic') {
    # Matrix where bootstraps will be stored
    pars.boot <- matrix(NA,bootstrap_iter,length(model$param_values))
    # Lower and Upper bounds
    L <- as.numeric(model$bounds["L", ])
    U <- as.numeric(model$bounds["U", ])
    cat("Starting bootstrap MCMSEM\n")
    pb <- txtProgressBar(0, bootstrap_iter, style = 3, width=min(c(options()$width, 107)))
      }
    if (se_type == 'one-step') {
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
      SEs <- apply(pars.boot, 2, sd)
      timeend <- Sys.time()
      boot1 <- timeend-timestart
    } else if (se_type == 'two-step') {
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
      SEs <- apply(pars.boot2, 2, sd)
      endtime <- Sys.time()
      boot2 <- endtime - starttime
    }
    #cat("\n")  # Prevent things from printing over completed progress bar
    # Table of point estimates and SE's
    results <- rbind(results, c(SEs))
  }

  colnames(results) <- model$param_names
  rownames(results) <- if(compute_se) c("est", "se") else c("est")
  return(results)
}

