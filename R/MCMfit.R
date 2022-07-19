# Fit MCM model
# Function to quadratically scale loss if estimates are out of bounds
.loss_power_bounds <- function(par, torch_bounds) {
  # Does not work for some reason, loss does not change, then jumps to NA
  lbound_check <- torch_less_equal(torch_cat(par), torch_bounds[["L"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  ubound_check <- torch_greater_equal(torch_cat(par), torch_bounds[["U"]])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  lbound_dist <- torch_square(torch_cat(par) - torch_bounds[["L"]])  # Distance from lbound
  ubound_dist <- torch_square(torch_cat(par) - torch_bounds[["U"]])  # Distance from ubound
  ldist_sum <- torch_sqrt(torch_sum(lbound_dist * lbound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  udist_sum <- torch_sqrt(torch_sum(ubound_dist * ubound_check)+1e-15)  # Distance multiplied by check: 0 if in bounds
  # 1e-15 is added to the lines above to prevent NAN gradients due to sqrt(0) issue
  pow <- torch_sum(ldist_sum) + torch_sum(udist_sum)
  return(pow)
}
# Objective function
.torch_objective <- function(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis, low_memory) {
  A <- torch_add(base_matrices[['A']], torch_sum(torch_mul(torch_maps[['A']], .par_list[['A']]), dim=3))
  # A <- torch_add(base_matrices[['A']], torch_sparse_coo_tensor(torch_maps[['A']]$indices()[1:2,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda')), .par_list[['A']][torch_maps[['A']]$indices()[3,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda'))], torch_maps[['A']]$shape[1:2]))
  Fm <- torch_add(base_matrices[['Fm']], torch_sum(torch_mul(torch_maps[['Fm']], .par_list[['Fm']]), dim=3))
  # Fm <- torch_add(base_matrices[['Fm']], torch_sparse_coo_tensor(torch_maps[['Fm']]$indices()[1:2,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda')), .par_list[['Fm']][torch_maps[['Fm']]$indices()[3,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda'))], torch_maps[['Fm']]$shape[1:2]))
  S <- torch_add(base_matrices[['S']], torch_sum(torch_mul(torch_maps[['S']], .par_list[['S']]), dim=3))
  # S <- torch_add(base_matrices[['S']], torch_sparse_coo_tensor(torch_maps[['S']]$indices()[1:2,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda')), .par_list[['S']][torch_maps[['S']]$indices()[3,] + torch_tensor(1, dtype=torch_long(), device=torch_device('cuda'))], torch_maps[['S']]$shape[1:2]))

  if (use_skewness) {Sk <- torch_add(base_matrices[['Sk']], torch_sum(torch_mul(torch_maps[['Sk']], .par_list[['Sk']]), dim=3))}
  if (use_kurtosis) {K <- torch_add(torch_add(torch_mul(torch_matmul(torch_matmul(torch_sqrt(S), base_matrices[['K']]), .torch_kron(torch_sqrt(S), .torch_kron(torch_sqrt(S), torch_sqrt(S)))), torch_masks[['K']]), torch_sum(torch_mul(torch_maps[['K']], .par_list[['K']]), dim=3)), base_matrices[['K2']])}
  #if (use_skewness) {Sk <- torch_add(base_matrices[['Sk']],torch_transpose(torch_mul(torch_hstack(list(torch_maps[['Sk_part1']], .par_list[['Sk']])), torch_maps[['Sk_part2']]), 1, 2))}
  #if (use_kurtosis) {K <- torch_add(torch_add(torch_mul(torch_matmul(torch_matmul(torch_sqrt(S), base_matrices[['K']]), .torch_kron(torch_sqrt(S), .torch_kron(torch_sqrt(S), torch_sqrt(S)))), torch_masks[['K']]), torch_transpose(torch_mul(torch_hstack(list(torch_maps[['K_part1']], .par_list[['K']])), torch_maps[['K_part2']]), 1, 2)), base_matrices[['K2']])}

  # Rstyle: M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), S), torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), torch_transpose(Fm, 1, 2))
  if (use_skewness) {
    # Rstyle: M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
    M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), Sk), .torch_kron(torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2))), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
  }
  if (use_kurtosis) {
    # Rstyle: M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))
    M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, torch_inverse(base_matrices[['diag_n_p']] - A)), K), torch_matmul(.torch_kron(.torch_kron( torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)),  torch_transpose(torch_inverse(base_matrices[['diag_n_p']] - A), 1, 2)), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))
  }
  ### Loss function
  # Rstyle: value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  if (use_skewness & use_kurtosis) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) +
      lossfunc(torch_mul(M3, m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']])) +
      lossfunc(torch_mul(M4, m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  } else if (use_skewness) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) +
      lossfunc(torch_mul(M3, m2v_masks[['m3']]), torch_mul(M3.obs, m2v_masks[['m3']]))
  } else if (use_kurtosis) {
    value <- lossfunc(torch_mul(M2, m2v_masks[['m2']]), torch_mul(M2.obs, m2v_masks[['m2']])) +
      lossfunc(torch_mul(M4, m2v_masks[['m4']]), torch_mul(M4.obs, m2v_masks[['m4']]))
  } else {
    stop("Oops, something went wrong.")
  }
  if (use_bounds) {
    pow <- .loss_power_bounds(.par_list, torch_bounds)
    loss <- value * (2^pow)
    return(loss)
  } else {
    return(value)
  }
}

# Fit wrapper function
.torch_fit <- function(M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, torch_dtype, return_history=FALSE, low_memory) {
  loss_hist <- c()
  lossfunc <- nn_mse_loss(reduction='sum')
  optim <- optim_rprop(.par_list,lr = learning_rate)
  if (low_memory) {gc(verbose=FALSE, full=TRUE)}
  for (i in 1:optim_iters[1]) {
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis)
    if (low_memory) {gc(verbose=FALSE, full=TRUE)}
    loss$backward()
    loss_hist <- c(loss_hist, loss$detach())
    optim$step()
    if (!(silent)) cat(paste0("loss", as.numeric(loss), "\n"))
  }
  calc_loss_torchfit <- function() {
    optim$zero_grad()
    loss <- .torch_objective(.par_list, lossfunc, torch_bounds, torch_masks, torch_maps, base_matrices, M2.obs, M3.obs, M4.obs, m2v_masks, use_bounds, use_skewness, use_kurtosis)
    if (low_memory) {gc(verbose=FALSE, full=TRUE)}
    loss$backward()
    loss_hist <<- c(loss_hist, loss$detach())
    if (!(silent)) {cat(paste0("loss", as.numeric(loss), "\n"))}
    return(loss)
  }
  # Use lbfgs to get really close....
  learning_rate <- 1
  optim <- optim_lbfgs(.par_list,lr= learning_rate)
  for (i in 1:optim_iters[2]) {
    optim$step(calc_loss_torchfit)
  }
  if (return_history) {
    return(list(par=.par_list, loss_hist=loss_hist))
  } else {
    return(.par_list)
  }

}

# Exported MCMfit function
MCMfit <- function(mcmmodel, data, compute_se=TRUE, se_type='asymptotic', optim_iters=c(50, 12), bootstrap_iter=200,bootstrap_chunks=1000,
                   learning_rate=0.02, silent=TRUE, use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE, device=NULL, low_memory=FALSE) {
  START_MCMfit <- Sys.time()
  model <- mcmmodel$copy()  # Model is changed if either use_skewness or use_kurtosis is set to FALSE, so I make a local copy here to ensure the original object stays intact
  if (is.null(device)) {
    device <- torch_device("cpu")
  }
  if (device == torch_device("cuda")) {
    if (!(silent)) {
      warning("silent was re-enabled as it is incompatible with a CUDA device.")
      silent <- TRUE
    }
  }
  cpu_device <- torch_device("cpu")
  if (!(se_type %in% c('two-step', 'one-step','asymptotic')))
    stop("se_type should be one of c('two-step', 'one-step','asymptotic')")
  if (ncol(data) != 2)
    warning("Use of a dataframe with more than 2 columns is still experimental.")
  if (!(use_skewness) & !(use_kurtosis))
    stop("At least either skewness or kurtosis has to be used")
  if (length(optim_iters) != 2) {
    if (length(optim_iters) == 2) {
      optim_iters <- c(optim_iters, 12) # If one value is passed to optim_iters I'm assuming that would be for RPROP
    } else {
      stop("optim_iters should be of length 2")
    }
  }
  # This code remains so this can be easily changed in future versions. As of torch 0.8.0 inverse does not work
  #  with half precision, so for now this is locked at 32 bit
  torch_precision <- 32
  if (torch_precision == 16) {
    torch_dtype <- torch_float16()
  } else if (torch_precision == 32) {
    torch_dtype <- torch_float32()
  } else {
    torch_dtype <- torch_float64()
  }
  if (nrow(data) < 1000)
    stop("Currently only a dataframe with at least 1000 rows is supported.")
  if (ncol(data) != model$meta_data$n_phenotypes) {
    msg <- paste0("Model expected ", model$meta_data$n_phenotypes, " phenotypes but ", ncol(data), " were found in the data. Please create a new model with this dataset.")
    stop(msg)
  }

  data <- as.matrix(data)
  # Scale data
  if ((model$meta_data$scale_data) & !(model$meta_data$data_was_scaled)) {
      data <- apply(data, 2, scale)
  }

  if (compute_se) ###????
    model_copy <- model$copy()
  if (!(use_kurtosis)) {
    kurt_par_idx <- which(startsWith(model$param_names, "k"))
    model$start_values <- model$start_values[, -kurt_par_idx]
    model$bounds <- model$bounds[, -kurt_par_idx]
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
    model$start_values <- model$start_values[, -skew_par_idx]
    model$bounds <- model$bounds[, -skew_par_idx]
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
  # Obtain covariance, coskewness and cokurtosis matrices
  M2.obs <- torch_tensor(cov(data), device=device, dtype=torch_dtype)
  M3.obs <- torch_tensor(M3.MM(data), device=device, dtype=torch_dtype)
  M4.obs <- torch_tensor(M4.MM(data), device=device, dtype=torch_dtype)

  torch_matrices <- .get_torch_matrices(model, device, M2.obs, M3.obs, M4.obs, torch_dtype)
  param_list <- torch_matrices[['param_list']]
  m2v_masks <- torch_matrices[['m2v_masks']]
  torch_bounds <- torch_matrices[['torch_bounds']]
  torch_masks <- torch_matrices[['torch_masks']]
  torch_maps <- torch_matrices[['torch_maps']]
  base_matrices <- torch_matrices[['base_matrices']]
  .par_list <- torch_matrices[['.par_list']]
  if (low_memory) {
    rm(torch_matrices)
    gc(full=TRUE)
  }
  # Obtain estimates with optimizer
  START_optim <- Sys.time()
  TIME_prep <-  START_optim - START_MCMfit
  out <- .torch_fit(M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, torch_dtype, return_history = TRUE, low_memory=low_memory)
  .par_tensor <- out[['par']]
  loss_hist <- as.numeric(torch_tensor(torch_vstack(out[["loss_hist"]]), device=cpu_device))
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
      for (i in 1:bootstrap_iter){
        ## Progress bar stuff
        setTxtProgressBar(pb, i)
        model2 <- model_copy$copy()  # Create new empty model
        #1. Sample from data with replacement
        boot <-   sample(1:nrow(data),nrow(data),T)
        sample <- data[boot,]

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
        .par_tensor <- .torch_fit(M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, torch_dtype, low_memory=low_memory)
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
        # 3. Sample cov/cosk/cokrt matrices and obtain mean of the sampled matrices to use as cov/cosk/cokrt matrix in the model
        boot <-   sample(1:bootstrap_chunks,bootstrap_chunks,T)
        M2.obs <-   torch_tensor(matrix(colMeans(sample.cov [boot,-1]),ncol(data),ncol(data), byrow=T), device=device, dtype=torch_dtype)
        M3.obs <-   torch_tensor(matrix(colMeans(sample.cosk[boot,-1]),ncol(data),ncol(data)^2, byrow=T), device=device, dtype=torch_dtype)
        M4.obs <-   torch_tensor(matrix(colMeans(sample.cokr[boot,-1]),ncol(data),ncol(data)^3, byrow=T), device=device, dtype=torch_dtype)
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
        .par_tensor <- .torch_fit(M2.obs, M3.obs, M4.obs, m2v_masks, torch_bounds, torch_masks, torch_maps, base_matrices, .par_list, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, torch_dtype, low_memory=low_memory)
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
  for (i in 1:length(model$param_coords)) {
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
  return(mcmresultclass(df=results, loss=loss_hist[length(loss_hist)], model=model$copy(),
                        history=list(
                          loss=loss_hist
                        ),
                        runtimes=list(Preparation=TIME_prep, Optimizer=TIME_optim, SE=TIME_se, Total=TIME_total),
                        info=list(version="0.6.1", compute_se=compute_se, se_type=se_type, optim_iters=optim_iters,
                                  bootstrap_iter=bootstrap_iter,bootstrap_chunks=bootstrap_chunks, learning_rate=learning_rate,
                                  silent=silent, use_bounds=use_bounds, use_skewness=use_skewness, use_kurtosis=use_kurtosis,
                                  device=device$type, low_memory=low_memory
                        )
                      ))
}

