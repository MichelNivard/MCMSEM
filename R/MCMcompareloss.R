MCMcompareloss <- function(results, data, weights=NULL, loss_type='auto', extensive_model_info=FALSE) {
  if (is.null(names(results))) {
    names(results) <- paste0("model", seq_len(length(results)))
  }
  results_ <- list()
  for (i in names(results)) {
    if (class(results[[i]])[[1]] == "mcmresultclass") {
      results_[[i]] <- results[[i]]
    }
  }
  results <- results_
  if (length(results) < 1) {
    stop("Results should be a list of at least one MCM result object")
  }
  if (class(data)[[1]] != "mcmdataclass") {
    if (!(class(data)[[1]] %in% c("data.frame", "matrix"))) {
      stop("Data should be an MCM data summary, dataframe, or matrix")
    }
    # If raw data is provided, check if input to weights argument matches that used in results
    for (i in names(results)) {
      if (!xor(is.null(weights), results[[i]]$info$weighted)) {
        if (is.null(weights)) {stop(paste0("No weights provided, but result ",i," was produced with weights."))} else {stop(paste0("Weights provided but result ",i," was produced without weights"))}
      }
    }
    data <- MCMdatasummary(data, scale_data=results[[1]]$model$meta_data$scale_data, weights=weights, prep_asymptotic_se=FALSE)
  } else {
    # If a data summary is provided, check if its weighted setting matches that used in results
    for (i in names(results)) {
      if (xor(data$meta_data$weighted, results[[i]]$info$weighted)) {
        if (!(data$meta_data$weighted)) {stop(paste0("Data is not weighted but result ",i," was produced with weights."))} else {stop(paste0("Data is weighted but result ",i," was produced without weights"))}
      }
    }
  }
  # Determine loss function(s) to use
  losses_used <- c()
  for (i in names(results)) {
    losses_used <- c(losses_used, results[[i]]$info$loss_type)
  }
  losses_used <- unique(losses_used)
  if (length(losses_used) > 1) {
    warning("Different loss functions were used for results provided.")
  }
  if (length(loss_type) == 1) {
    if (loss_type == 'auto') {
      losses <- losses_used
    } else if (!(loss_type %in% c("mse", "smooth_l1" ,'l1'))) {
      stop("loss_type should be one of c('auto', 'mse', 'smooth_l1') or a vector with a combination of these")
    } else {
      losses <- c(loss_type)
    }
  } else {
    losses <- c()
    for (i in loss_type) {
      if (!(i %in% c("mse", "smooth_l1" ,'l1'))) {
        cat(paste0("loss type '",i,"' in vector not recognized, and therefore ignored.\n"))
      } else {
        losses <- c(losses, i)
      }
    }
    losses <- unique(losses)
  }
  out <- data.frame(row.names=names(results))

  lossfuncs <- list(
    mse=nn_mse_loss(reduction='sum'),
    smooth_l1=nn_smooth_l1_loss(reduction='sum'),
    l1=nn_l1_loss(reduction='sum')
  )
  m2v_masks <- list(
    m2=.torch_m2m2v_mask(torch_tensor(data$M2), device=torch_device("cpu"), dtype=torch_float32()),
    m3=.torch_m3m2v_mask(torch_tensor(data$M3), device=torch_device("cpu"), dtype=torch_float32()),
    m4=.torch_m4m2v_mask(torch_tensor(data$M4), device=torch_device("cpu"), dtype=torch_float32())
  )
  # Check if models were fitted on the same data
  if (length(results) > 1) {
    M2_model1 <- results[[names(results)[1]]]$history$M2obs
    for (resname in names(results)[2:length(results)]) {
      if (sum((round(results[[resname]]$history$M2obs - M2_model1, 3)))!=0) {
        warning("Results seem to have been create using different datasets")
      }
    }
  }
  # Determine loss for each model, and difference with model1 loss for models 2-N
  comoments_used <- c()
  for (loss_type in losses) {
    lossfunc <- lossfuncs[[loss_type]]
    i <- 0
    for (resname in names(results)) {
      if (results[[resname]]$info$use_kurtosis & results[[resname]]$info$use_skewness) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$history$M2pred, M3=results[[resname]]$history$M3pred, M4=results[[resname]]$history$M4pred), m2v_masks, data$M2, data$M3, data$M4, TRUE, TRUE)
        comoments_used <- c(comoments_used, "kurtskewvar")
      } else if (results[[resname]]$info$use_skewness) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$history$M2pred, M3=results[[resname]]$history$M3pred), m2v_masks, data$M2, data$M3, data$M4, TRUE, FALSE)
        comoments_used <- c(comoments_used, "skewvar")
      } else if (results[[resname]]$info$use_kurtosis) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$history$M2pred, M4=results[[resname]]$history$M4pred), m2v_masks, data$M2, data$M3, data$M4, FALSE, TRUE)
        comoments_used <- c(comoments_used, "kurtvar")
      } else {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$history$M2pred), m2v_masks, data$M2, data$M3, data$M4, FALSE, FALSE)
        comoments_used <- c(comoments_used, "var")
      }
      loss <- as.numeric(loss)
      out[resname, paste0(loss_type,"_loss")] <- loss
      if (i > 0) {
        out[resname, paste0(loss_type,"_diff")] <- out[1, paste0(loss_type,"_loss")] - loss
      }
      i <- i + 1
    }
  }
  if (length(unique(comoments_used)) > 1) {
    stop("Different co-moment matrices are used across models.")
  }
  # Store chisq and BIC for each model
  for (loss_type in losses) {
    for (resname in names(results)) {
      out[resname, paste0(loss_type,"_chisq")] <- data$meta_data$N * out[resname, paste0(loss_type,"_loss")]
      out[resname, paste0(loss_type,"_bic")] <- data$meta_data$N * out[resname, paste0(loss_type,"_loss")] + ncol(results[[resname]]$df) * log(data$meta_data$N)
    }
  }
  # Store N parameters and optim iters and learning rates for each model
  for (resname in names(results)) {
    out[resname, "N_parameters"] <- ncol(results[[resname]]$df)
    out[resname, "optim1_iters"] <- results[[resname]]$info$optim_iters[1]
    out[resname, "optim1_lr"] <- results[[resname]]$info$learning_rate[1]
    out[resname, "optim2_iters"] <- results[[resname]]$info$optim_iters[2]
    out[resname, "optim2_lr"] <- results[[resname]]$info$learning_rate[2]
  }
  # Store loss type, learning rate, use skew, use kurt
  for (resname in names(results)) {
    for (i in c("loss_type", "use_skewness", "use_kurtosis")) {
      out[resname, i] <- results[[resname]]$info[[i]]
    }
  }
  if (!extensive_model_info) {
    # If all optim_iters, learning rates, use skew, use kurt are equal they are not informative so can be removed
    if ((length(unique(out[, "optim1_iters"]))==1) & (length(unique(out[, "optim2_iters"]))==1)) {out[,c("optim1_iters", "optim2_iters")] <- NULL}
    if ((length(unique(out[, "optim1_lr"]))==1) & (length(unique(out[, "optim2_lr"]))==1)) {out[,c("optim1_lr", "optim2_lr")] <- NULL}
    for (i in c("loss_type", "use_skewness", "use_kurtosis")) {
      if (length(unique(out[, i]))==1) {out[,i] <- NULL}
    }
  } else {
    # Store version, use_bounds, optimizer, compute se, se type, silent, device and low memory
    for (resname in names(results)) {
      for (i in c("version", "use_bounds",  "optimizer", "compute_se", "se_type", "silent", "device", "low_memory")) {
        out[resname, i] <- results[[resname]]$info[[i]]
      }
    }
    # Store runtimes
    for (resname in names(results)) {
      for (i in c("Preparation", "Optimizer", "SE", "Total")) {
        out[resname, paste0(i,"_runtime")] <- as.numeric(results[[resname]]$runtimes[[i]], units='secs')
      }
    }
  }
  return(out)
}