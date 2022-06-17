# Fit MCM model
# Objective function
# torch-ready versions of
.loss_power_bounds <- function(.par, model) {
  par <- as.numeric(.par)
  bounds <- as.matrix(model$bounds)
  lbound_check <- par > bounds[1,]
  ubound_check <- par < bounds[2,]
  if (all(lbound_check & ubound_check)) {
    return(0)
  } else {
    pow <- 1
    if (!all(lbound_check)) {
      idx_fail <- which(!(lbound_check))
      # Add absolute distance from the bounds to power
      pow <- pow + sum(abs(par[idx_fail]) - abs(bounds[1, idx_fail]))
    }
    if (!all(ubound_check)) {
      idx_fail <- which(!(lbound_check))
      # Add absolute distance from the bounds to power
      pow <- pow + sum(abs(par[idx_fail]) - abs(bounds[1, idx_fail]))
    }
    return(pow)
  }
}

.loss_power_bounds_torch <- function(.par, model) {
  # Does not work for some reason, loss does not change, then jumps to NA
  torch_bounds <- torch_tensor(as.matrix(model$bounds))  # If it works, this only needs to be done once and can be passed as argument
  lbound_check <- !(.par > torch_bounds[1,])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  ubound_check <- !(.par < torch_bounds[2,])  # 0 if parameter is in bounds, 1 if parameter is out of bounds
  lbound_dist <- torch_square(.par - torch_bounds[1, ])  # Distance from lbound
  ubound_dist <- torch_square(.par - torch_bounds[2, ])  # Distance from ubound
  ldist_sum <- torch_sqrt(torch_sum(lbound_dist * lbound_check))  # Distance multiplied by check: 0 if in bounds
  udist_sum <- torch_sqrt(torch_sum(ubound_dist * ubound_check))  # Distance multiplied by check: 0 if in bounds
  return(torch_sum(ldist_sum) + torch_sum(udist_sum))
}

.torch_objective <- function(.par, model, torch_coords, M2.obs, M3.obs, M4.obs, use_bounds, use_skewness, use_kurtosis) {
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
  if (use_kurtosis) {
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
  }

  ###### Compute the observed cov, cosk, and cokurt matrices #################
  diag_np_a_inv <- torch_inverse(diag_n_p - A)
  diag_np_a_inv_t <- torch_transpose(torch_inverse(diag_n_p - A), 1, 2)
  # M2 <- Fm %*% solve(diag(n_p) - A) %*% S %*%  t(solve(diag(n_p)-A))  %*% t(Fm)
  M2 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), S), diag_np_a_inv_t), torch_transpose(Fm, 1, 2))
  if (use_skewness) {
    # M3 <- Fm %*% solve(diag(n_p) - A) %*% Sk %*% (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm))
    M3 <- torch_matmul(torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), Sk), .torch_kron(diag_np_a_inv_t, diag_np_a_inv_t)), .torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)))
  }
  if (use_kurtosis) {
    # M4 <- Fm %*% solve(diag(n_p) - A) %*% K %*%  (t(solve(diag(n_p)-A)) %x% t(solve(diag(n_p)-A))  %x%  t(solve(diag(n_p)-A))) %*% (t(Fm) %x% t(Fm) %x% t(Fm))
    M4 <- torch_matmul(torch_matmul(torch_matmul(Fm, diag_np_a_inv), K), torch_matmul(.torch_kron(.torch_kron(diag_np_a_inv_t, diag_np_a_inv_t), diag_np_a_inv_t), .torch_kron(.torch_kron(torch_transpose(Fm, 1, 2), torch_transpose(Fm, 1, 2)), torch_transpose(Fm, 1, 2))))
  }

  ### Loss function
  # value <- sum((.m2m2v(M2.obs) - .m2m2v(M2))^2) + sum((.m3m2v(M3.obs) - .m3m2v(M3))^2) + sum((.m4m2v(M4.obs)- .m4m2v(M4))^2)
  # value <- torch_sum((torch_square(torch_tril(M2.obs) - torch_tril(M2))))  +  torch_sum((torch_square(torch_tril(M3.obs$reshape(c(2,2,2))) - torch_tril(M3$reshape(c(2,2,2)))))) + torch_sum((torch_square(torch_tril(M4.obs$reshape(c(2,2,2,2)))- torch_tril(M4$reshape(c(2,2,2,2))))))
  if (use_skewness & use_kurtosis) {
    value <- torch_sum((torch_square(.torch_m2m2v(M2.obs) - .torch_m2m2v(M2))))  +  torch_sum((torch_square(.torch_m3m2v(M3.obs) - .torch_m3m2v(M3)))) + torch_sum((torch_square(.torch_m4m2v(M4.obs)- .torch_m4m2v(M4))))
  } else if (use_skewness) {
    value <- torch_sum((torch_square(.torch_m2m2v(M2.obs) - .torch_m2m2v(M2))))  +  torch_sum((torch_square(.torch_m3m2v(M3.obs) - .torch_m3m2v(M3))))
  } else if (use_kurtosis) {
    value <- torch_sum((torch_square(.torch_m2m2v(M2.obs) - .torch_m2m2v(M2)))) + torch_sum((torch_square(.torch_m4m2v(M4.obs)- .torch_m4m2v(M4))))
  } else {
    stop("Oops, something went wrong.")
  }
  if (use_bounds) {
    pow <- .loss_power_bounds(.par, model)
    loss <- value * (2^pow)
    return(loss)
  } else {
    return(value)
  }

}



.torch_fit <- function(model, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, return_history=FALSE) {
  .par <- torch_tensor(as.numeric(model$start_values), requires_grad=TRUE)
  loss_hist <- c()
  optim <- optim_rprop(.par,lr = learning_rate)
  for (i in 1:optim_iters[1]) {
    optim$zero_grad()
    loss <- .torch_objective(.par, model, torch_coords, M2.obs, M3.obs, M4.obs, use_bounds, use_skewness, use_kurtosis)
    loss_hist <- c(loss_hist, as.numeric(loss))
    loss$backward()
    optim$step()
    if (!(silent)) cat(paste0("loss", as.numeric(loss), "\n"))
  }
  calc_loss_torchfit <- function() {
    optim$zero_grad()
    loss <- .torch_objective(.par, model, torch_coords, M2.obs, M3.obs, M4.obs, use_bounds, use_skewness, use_kurtosis)
    loss_hist <<- c(loss_hist, as.numeric(loss))
    if (!(silent)) {cat(paste0("loss", as.numeric(loss), "\n"))}
    loss$backward()
    return(loss)
  }
  # Use lbfgs to get really close....
  learning_rate <-1
  optim <- optim_lbfgs(.par,lr= learning_rate)
  for (i in 1:optim_iters[2]) {
    optim$step(calc_loss_torchfit)
  }
  if (return_history) {
    return(list(par=.par, loss_hist=loss_hist))
  } else {
    return(.par)
  }

}

MCMfit <- function(mcmmodel, data, compute_se=TRUE, se_type='asymptotic', optim_iters=c(50, 12), bootstrap_iter=200,bootstrap_chunks=1000,
                   learning_rate=0.02, silent=TRUE, use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE) {
  model <- mcmmodel$copy()  # Model is changed if either use_skewness or use_kurtosis is set to FALSE, so I make a local copy here to ensure the original object stays intact
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
  M2.obs <- torch_tensor(cov(data))
  M3.obs <- torch_tensor(M3.MM(data))
  M4.obs <- torch_tensor(M4.MM(data))

  torch_coords <- .get_torch_coords(model)
  # Obtain estimates with optimizer
  out <- .torch_fit(model, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis, return_history = TRUE)
  .par <- out[['par']]
  loss_hist <- out[['loss_hist']]
  # Store estimates including minimization objective, using this to evaluate/compare fit
  results <-  as.data.frame(matrix(as.numeric(.par), nrow = 1))

  if (compute_se) {

    if(se_type == 'asymptotic'){
      SEs <- .std.err(data=data,par=as.numeric(.par),model=model, use_skewness, use_kurtosis)
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
        M2.obs <- torch_tensor(cov(sample))
        M3.obs <- torch_tensor(M3.MM(sample))
        M4.obs <- torch_tensor(M4.MM(sample))

        #3. Fit model
        # Estimate parameters with model function specified above
        .par <- .torch_fit(model2, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis)
        # Store point estimates of bootstraps
        pars.boot[i,] <- as.numeric(.par)

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
        M2.obs <-   torch_tensor(matrix(colMeans(sample.cov [boot,-1]),ncol(data),ncol(data), byrow=T))
        M3.obs <-   torch_tensor(matrix(colMeans(sample.cosk[boot,-1]),ncol(data),ncol(data)^2, byrow=T))
        M4.obs <-   torch_tensor(matrix(colMeans(sample.cokr[boot,-1]),ncol(data),ncol(data)^3, byrow=T))

        #4. Fit model
        # Estimate parameters with model function specified above
        .par <- .torch_fit(model2, torch_coords, M2.obs, M3.obs, M4.obs, learning_rate, optim_iters, silent, use_bounds, use_skewness, use_kurtosis)
        # Store point estimates of bootstraps
        pars.boot[i,] <- as.numeric(.par)

      }
      close(pb)
      SEs <- apply(pars.boot2, 2, sd)
    }
    # Table of point estimates and SE's
    results <- rbind(results, c(SEs))
  }

  colnames(results) <- model$param_names
  rownames(results) <- if(compute_se) c("est", "se") else c("est")
  return(mcmresultclass(df=results, loss=loss_hist[length(loss_hist)],
                      history=list(
                        model=model$copy(),
                        loss=loss_hist
                      )))
}

