# Fit MCM model
# Local functions used in MCMfit moved to local_fit.R

# Exported MCMfit function
MCMfit <- function(mcmmodel, data, weights=NULL, compute_se=TRUE, se_type='asymptotic', optimizer="rprop", optim_iters=c(50, 12), loss_type='mse', bootstrap_iter=200,bootstrap_chunks=1000,
                   learning_rate=c(0.02, 1), silent=TRUE, use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE, device=NULL, low_memory=FALSE, outofbounds_penalty=1, debug=FALSE) {
  START_MCMfit <- Sys.time()
  if (debug) {cat("Verifying input...\n")}
  if (class(mcmmodel)[[1]] == "mcmresultclass") {
    cat("Note: Using the old model from the result object provided\n")
    # If a result class is provided, first check if settings are similar enough
    if (use_kurtosis & !(mcmmodel$info$use_kurtosis)) {cat("      use_kurtosis reset to FALSE for the current fit, as it was set to FALSE previously. \n"); use_kurtosis <- FALSE}
    if (use_skewness & !(mcmmodel$info$use_skewness)) {cat("      use_skewness reset to FALSE for the current fit, as it was set to FALSE previously. \n"); use_skewness <- FALSE}
    if (class(data)[[1]] != "mcmdataclass") {cat("      Next time when you intend to fit your model twice, make an MCMdata summary first. That will save you time. \n")}
    model <- mcmmodel$model$copy()
    model$start_values$set_all(mcmmodel$model$param_values)
  } else {
    model <- mcmmodel$copy()  # Model is changed if either use_skewness or use_kurtosis is set to FALSE, so I make a local copy here to ensure the original object stays intact
    if (!(use_kurtosis)) {
      kurt_par_idx <- which(startsWith(model$param_names, "k"))
      model$start_values$drop(kurt_par_idx)
      model$param_names <- model$param_names[-kurt_par_idx]
      model$param_values <- model$param_values[-kurt_par_idx]
      new_coords <- list()
      iter <- 1
      for (i in model$param_coords) {
        if (!(i[[1]] == "K")) {
          new_coords[[iter]] <- i
          iter <- iter +1
        }
      }
      model$param_coords <- new_coords
    } else if (!(use_skewness)) {
      skew_par_idx <- which(startsWith(model$param_names, "sk"))
      model$start_values$drop(skew_par_idx)
      model$param_names <- model$param_names[-skew_par_idx]
      model$param_values <- model$param_values[-skew_par_idx]
      new_coords <- list()
      iter <- 1
      for (i in model$param_coords) {
        if (!(i[[1]] == "Sk")) {
          new_coords[[iter]] <- i
          iter <- iter +1
        }
      }
      model$param_coords <- new_coords
    }
  }
  if (class(data)[[1]] != "mcmdataclass") {
    data_org <- data
    if (debug) {cat("Converting data\n")}
    data <- MCMdatasummary(data, scale_data=model$meta_data$scale_data, weights=weights, prep_asymptotic_se=((compute_se) & (se_type == "asymptotic")))
  } else {
    if (compute_se) {
      if (se_type %in% c("one-step", 'two-step'))
        stop("Bootstrap SE computation not supported with summary data")
      if ((se_type == 'asymptotic') & !(data$SE$computed))
        stop("This summary data was made without the required preparation for asymptotic SE")
    }
  }
  model$meta_data$n_obs <- data$meta_data$n # Store the N of the training data, note this may not always be the same as the N the model was initially made with
  if (is.null(device)) {
    device <- torch_device("cpu")
  }
  if (device == torch_device("cuda")) {
    if (!(silent)) {
      warning("silent was re-enabled as it is incompatible with a CUDA device.")
      silent <- TRUE
    }
  }
  if (model$meta_data$weighted & is.null(weights)) {stop("Model was created with weights but no weights are provided")}
  cpu_device <- torch_device("cpu")
  if (!(se_type %in% c('two-step', 'one-step','asymptotic')))
    stop("se_type should be one of c('two-step', 'one-step','asymptotic')")
  if (data$meta_data$ncol != 2)
    warning("Use of a dataframe with more than 2 columns is still experimental.")
  if (!(use_skewness) & !(use_kurtosis))
    stop("At least either skewness or kurtosis has to be used")
  if (length(optim_iters) != 2) {
    if (length(optim_iters) == 1) {
      optim_iters <- c(optim_iters, 12) # If one value is passed to optim_iters I'm assuming that would be for RPROP
    } else {
      stop("optim_iters should be of length 2")
    }
  }
  if (length(learning_rate) != 2) {
    if (length(learning_rate) == 1) {
      learning_rate <- c(learning_rate, 1.0) # If one value is passed to optim_iters I'm assuming that would be for RPROP
    } else {
      stop("learning_rate should be of length 2")
    }
  }
  if (outofbounds_penalty < 0) {
    stop("outofbounds_penalty should be a positive integer")
  } else if ((outofbounds_penalty == 0) & (use_bounds)) {
    cat("Setting use_bounds to FALSE as this is more efficient and has the same effect as outofbounds_penalty=0\n")
    use_bounds <- FALSE
  }
  # This code remains so this can be easily changed in future versions. As of torch 0.8.0 inverse does not work
  #  with half precision, so for now this is locked at 32 bit
  torch_precision <- 32
  if (torch_precision == 16) {torch_dtype <- torch_float16()} else if (torch_precision == 32) {torch_dtype <- torch_float32()} else if (torch_precision == 64) {torch_dtype <- torch_float64()} else {stop("Precision not recognized, should be one of (16, 32, 64)")}
  if (data$meta_data$ncol != model$meta_data$n_phenotypes) {
    msg <- paste0("Model expected ", model$meta_data$n_phenotypes, " phenotypes but ", data$meta_data$ncol, " were found in the data. Please create a new model with this dataset.")
    stop(msg)
  }

  lossfunc <- .get_lossfunc(loss_type)
  optimfunc <- .get_optimfunc(optimizer)

  if (compute_se)
    model_copy <- model$copy()
  if (debug) {cat("Transferring co-moments\n")}
  # Obtain covariance, coskewness and cokurtosis matrices
  M2.obs <- torch_tensor(data$M2, device=device, dtype=torch_dtype)
  M3.obs <- torch_tensor(data$M3, device=device, dtype=torch_dtype)
  M4.obs <- torch_tensor(data$M4, device=device, dtype=torch_dtype)
  if (debug) {cat("Generating torch matrices\n")}
  torch_matrices <- .get_torch_matrices(model, device, M2.obs, M3.obs, M4.obs, torch_dtype)
  param_list <- torch_matrices[['param_list']]
  m2v_masks <- torch_matrices[['m2v_masks']]
  torch_bounds <- torch_matrices[['torch_bounds']]
  torch_masks <- torch_matrices[['torch_masks']]
  torch_maps <- torch_matrices[['torch_maps']]
  base_matrices <- torch_matrices[['base_matrices']]
  .par_list <- torch_matrices[['.par_list']]
  if (low_memory) {
    if (debug) {cat("Clearing memory\n")}
    rm(torch_matrices)
    gc(full=TRUE)
  }
  # Obtain estimates with optimizer
  START_optim <- Sys.time()
  TIME_prep <-  START_optim - START_MCMfit
  if (debug) {cat("Strting fit\n")}
  out <- .torch_fit(optimfunc, M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, lossfunc, return_history = TRUE, low_memory=low_memory, outofbounds_penalty=outofbounds_penalty,debug=debug)
  .par_tensor <- out[['par']]
  loss_hist <- as.numeric(torch_tensor(torch_vstack(out[["loss_hist"]]), device=cpu_device))
  pred_matrices <- out[['pred_matrices']]
  # Store estimates including minimization objective, using this to evaluate/compare fit
  .par <- list()
  for (i in names(.par_tensor)) {
    if (length(param_list[[i]]) > 0) {
      .par[[i]] <- .par_tensor[[i]]  # Only include parameters that were actually estimated (requires_grad=TRUE)
    }
  }
  results <-  as.data.frame(matrix(as.numeric(torch_tensor(torch_cat(.par), device=cpu_device)), nrow = 1))
  START_se <- Sys.time()
  TIME_optim <- START_se - START_optim
  if (compute_se) {
    if(se_type == 'asymptotic'){
      # SEs <- .std.err(data=data,par=as.numeric(torch_tensor(torch_cat(.par), device=cpu_device)), model=model, use_skewness, use_kurtosis)
      SEs <- .std.err(data, .par_list, use_skewness, use_kurtosis, torch_masks, torch_maps, base_matrices, m2v_masks, device, low_memory)
    }

    if(se_type != 'asymptotic') {
    # Matrix where bootstraps will be stored
    pars.boot <- matrix(NA,bootstrap_iter,length(model$param_values))
    # Lower and Upper bounds
    cat("Starting bootstrap MCMSEM\n")
    pb <- txtProgressBar(0, bootstrap_iter, style = 3, width=min(c(options()$width, 107)))
      }
    if (se_type == 'one-step') {
      #### BOOT 1:  NORMAL BOOTSTRAP
      ##############################
      # Matrix where bootstraps will be stored
      # Bootstrap
      for (i in seq_len(bootstrap_iter)) {
        ## Progress bar stuff
        setTxtProgressBar(pb, i)
        model2 <- model_copy$copy()  # Create new empty model
        #1. Sample from data with replacement
        boot <-   sample(seq_len(nrow(data_org)), nrow(data_org), T)
        sample <- data_org[boot,]

        #2. Get covariance, coskewness and cokurtosis matrices
        M2.obs <- torch_tensor(cov(sample), device=device, dtype=torch_dtype)
        M3.obs <- torch_tensor(M3.MM(sample), device=device, dtype=torch_dtype)
        M4.obs <- torch_tensor(M4.MM(sample), device=device, dtype=torch_dtype)
        # Get new torch matrices from model copy
        torch_matrices <- .get_torch_matrices(model2, device, M2.obs, M3.obs, M4.obs, torch_dtype)
        m2v_masks <- torch_matrices[['m2v_masks']]
        torch_bounds <- torch_matrices[['torch_bounds']]
        torch_masks <- torch_matrices[['torch_masks']]
        torch_maps <- torch_matrices[['torch_maps']]
        base_matrices <- torch_matrices[['base_matrices']]
        .par_list <- torch_matrices[['.par_list']]
        # Obtain estimates with optimizer
        START_optim <- Sys.time()
        TIME_prep <-  START_optim - START_MCMfit

        #3. Fit model
        # Estimate parameters with model function specified above
        .par_tensor <- .torch_fit(optimfunc, M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, lossfunc, low_memory=low_memory, outofbounds_penalty=outofbounds_penalty)
        .par <- list()
        for (j in names(.par_tensor)) {
          if (length(param_list[[j]]) > 0) {
            .par[[j]] <- .par_tensor[[j]]  # Only include parameters that were actually estimated (requires_grad=TRUE)
          }
        }
        # Store point estimates of bootstraps
        pars.boot[i,] <- as.numeric(torch_tensor(torch_cat(.par), device=cpu_device))

      }
      close(pb)
      SEs <- apply(pars.boot, 2, sd)

    } else if (se_type == 'two-step') {
      #### BOOT 2: TWO STEP BOOTSTRAP
      ##############################
      ### STEP 1
      # 1. Bind the data to a random sample binning people into "bootstrap_chunks" groups
      step1 <- cbind(data,sample(1:bootstrap_chunks,nrow(data_org),replace=T))
      colnames(step1)[ncol(step1)] <- "group"
      step1 <- as.data.frame(step1)

      # 2. Get covariance, coskenwess and cokurtosis matrices per group
      sample.cov  <- aggregate(seq_len(nrow(step1)), by=list(step1$group), function(s) cov(step1[s, colnames(data_org)]))
      sample.cosk <- aggregate(seq_len(nrow(step1)), by=list(step1$group), function(s) M3.MM(as.matrix(step1[s, colnames(data_org)])))
      sample.cokr <- aggregate(seq_len(nrow(step1)), by=list(step1$group), function(s) M4.MM(as.matrix(step1[s, colnames(data_org)])))
      pars.boot2 <- matrix(NA,nrow=bootstrap_iter,ncol=length(model$param_values))

      ### STEP 2
      for (i in seq_len(bootstrap_iter)) {
        ## Progress bar stuff
        setTxtProgressBar(pb, i)
        model2 <- model_copy$copy()  # Create new empty model
        # 3. Sample cov/cosk/cokrt matrices and obtain mean of the sampled matrices to use as cov/cosk/cokrt matrix in the model
        boot <-   sample(1:bootstrap_chunks,bootstrap_chunks,T)
        M2.obs <-   torch_tensor(matrix(colMeans(sample.cov [boot,-1]),ncol(data_org),ncol(data_org), byrow=T), device=device, dtype=torch_dtype)
        M3.obs <-   torch_tensor(matrix(colMeans(sample.cosk[boot,-1]),ncol(data_org),ncol(data_org)^2, byrow=T), device=device, dtype=torch_dtype)
        M4.obs <-   torch_tensor(matrix(colMeans(sample.cokr[boot,-1]),ncol(data_org),ncol(data_org)^3, byrow=T), device=device, dtype=torch_dtype)
        # Get new torch matrices from model copy
        torch_matrices <- .get_torch_matrices(model2, device, M2.obs, M3.obs, M4.obs, torch_dtype)
        m2v_masks <- torch_matrices[['m2v_masks']]
        torch_bounds <- torch_matrices[['torch_bounds']]
        torch_masks <- torch_matrices[['torch_masks']]
        torch_maps <- torch_matrices[['torch_maps']]
        base_matrices <- torch_matrices[['base_matrices']]
        .par_list <- torch_matrices[['.par_list']]
        # Obtain estimates with optimizer
        START_optim <- Sys.time()
        TIME_prep <-  START_optim - START_MCMfit

        #3. Fit model
        # Estimate parameters with model function specified above
        .par_tensor <- .torch_fit(optimfunc, M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, lossfunc, low_memory=low_memory, outofbounds_penalty=outofbounds_penalty)
        .par <- list()
        for (j in names(.par_tensor)) {
          if (length(param_list[[j]]) > 0) {
            .par[[j]] <- .par_tensor[[j]]  # Only include parameters that were actually estimated (requires_grad=TRUE)
          }
        }
        # Store point estimates of bootstraps
        pars.boot2[i,] <- as.numeric(torch_tensor(torch_cat(.par), device=cpu_device))
      }
      close(pb)
      SEs <- apply(pars.boot2, 2, sd)
    }
    # Table of point estimates and SE's
    results <- rbind(results, SEs)
  }
  # Place resulting parameter estimates back into model matrices
  for (i in seq_along(model$param_coords)) {
    idat <- model$param_coords[[i]]
    model$num_matrices[[idat[[1]]]][idat[[2]]] <- as.numeric(results[1, i])
    # Currently, if there is a parameter like "-a1", this will still return the base a1 value.
    # If we want to change that so that the returned estimate is also flipped, use this
    # model$num_matrices[[idat[[1]]]][idat[[2]]] <- as.numeric(results[1, i]) * idat[[3]]
  }

  model$param_values <- as.numeric(results[1, ])
  colnames(results) <- model$param_names
  rownames(results) <- if(compute_se) c("est", "se") else "est"
  STOP <- Sys.time()
  TIME_se <- STOP - START_se
  TIME_total <- STOP - START_MCMfit
  history <- list(loss=loss_hist, M2pred=as.matrix(torch_tensor(pred_matrices$M2, device=cpu_device)), M2obs=as.matrix(torch_tensor(M2.obs, device=cpu_device)))
  if (use_skewness) {history[['M3pred']] <- as.matrix(torch_tensor(pred_matrices$M3, device=cpu_device)); history[['M3obs']] <- as.matrix(torch_tensor(M3.obs, device=cpu_device))}
  if (use_skewness) {history[['M4pred']] <- as.matrix(torch_tensor(pred_matrices$M4, device=cpu_device)); history[['M4obs']] <- as.matrix(torch_tensor(M4.obs, device=cpu_device))}
  loss <- .calc_loss(lossfunc, pred_matrices, m2v_masks, M2.obs, M3.obs, M4.obs, use_skewness, use_kurtosis)
                             
                               
  return(mcmresultclass(df=results, loss=as.numeric(torch_tensor(loss, device=cpu_device)), model=model$copy(),
                        history=history,
                        runtimes=list(Preparation=TIME_prep, Optimizer=TIME_optim, SE=TIME_se, Total=TIME_total),
                        info=list(version=MCMSEMversion, compute_se=compute_se, se_type=se_type, optim_iters=optim_iters,
                                  bootstrap_iter=bootstrap_iter,bootstrap_chunks=bootstrap_chunks, learning_rate=learning_rate,
                                  silent=silent, use_bounds=use_bounds, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                                  device=device$type, low_memory=low_memory, weighted=data$meta_data$weighted, loss_type=loss_type,
                                  optimizer=optimizer, n=data$meta_data$n)
                      ))
}

