MCMcompareloss <- function(results, test_data, weights=NULL, loss_type='auto', extensive_model_info=FALSE) {
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
  if (class(test_data)[[1]] != "mcmdataclass") {
    if (!(class(test_data)[[1]] %in% c("data.frame", "matrix"))) {
      stop("Data should be an MCM data summary, dataframe, or matrix")
    }
    # If raw data is provided, check if input to weights argument matches that used in results
    for (i in names(results)) {
      if (!xor(is.null(weights), results[[i]]$info$weighted)) {
        if (is.null(weights)) {stop(paste0("No weights provided, but result ",i," was produced with weights."))} else {stop(paste0("Weights provided but result ",i," was produced without weights"))}
      }
    }
    data <- MCMdatasummary(test_data, scale_data=results[[1]]$model$meta_data$scale_data, weights=weights, prep_asymptotic_se=FALSE)
  } else {
    data <- test_data
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
    } else {
      lossfunc_ <- .get_lossfunc(loss_type)  # This is just to produce an error if loss type is not recognizerd
      losses <- c(loss_type)
    }
  } else {
    losses <- c()
    for (i in loss_type) {
      lossfunc_ <- .get_lossfunc(i)  # This is just to produce an error if loss type is not recognizerd
      losses <- c(losses, i)
    }
    losses <- unique(losses)
  }
  out <- data.frame(row.names=names(results))

  m2v_masks <- list(
    m2=.torch_m2m2v_mask(torch_tensor(data$M2), device=torch_device("cpu"), dtype=torch_float32()),
    m3=.torch_m3m2v_mask(torch_tensor(data$M3), device=torch_device("cpu"), dtype=torch_float32()),
    m4=.torch_m4m2v_mask(torch_tensor(data$M4), device=torch_device("cpu"), dtype=torch_float32())
  )
  # Check if models were fitted on the same data
  if (length(results) > 1) {
    M2_model1 <- results[[names(results)[1]]]$observed$M2
    for (resname in names(results)[2:length(results)]) {
      if (sum((round(results[[resname]]$observed$M2 - M2_model1, 3)))!=0) {
        warning("Results seem to have been create using different datasets")
      }
    }
  }
  # Store original loss, chisq and bic for each model
  for (resname in names(results)) {
    out[resname, paste0(results[[resname]]$info$loss_type, "_train_loss")] <- results[[resname]]$loss
    out[resname, "train_chisq"] <- results[[resname]]$model$meta_data$n_obs * results[[resname]]$loss
    out[resname, "train_bic"] <- results[[resname]]$model$meta_data$n_obs * results[[resname]]$loss + ncol(results[[resname]]$df) * log(results[[resname]]$model$meta_data$n_obs)
  }

  # Determine loss for each model, and difference with model1 loss for models 2-N
  comoments_used <- c()
  for (loss_type in losses) {
    lossfunc <- .get_lossfunc(loss_type)
    i <- 0
    for (resname in names(results)) {
      if (results[[resname]]$info$use_kurtosis & results[[resname]]$info$use_skewness) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$predicted$M2, M3=results[[resname]]$predicted$M3, M4=results[[resname]]$predicted$M4), m2v_masks, data$M2, data$M3, data$M4, TRUE, TRUE)
        comoments_used <- c(comoments_used, "kurtskewvar")
      } else if (results[[resname]]$info$use_skewness) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$predicted$M2, M3=results[[resname]]$predicted$M3), m2v_masks, data$M2, data$M3, data$M4, TRUE, FALSE)
        comoments_used <- c(comoments_used, "skewvar")
      } else if (results[[resname]]$info$use_kurtosis) {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$predicted$M2, M4=results[[resname]]$predicted$M4), m2v_masks, data$M2, data$M3, data$M4, FALSE, TRUE)
        comoments_used <- c(comoments_used, "kurtvar")
      } else {
        loss <- .calc_loss(lossfunc, list(M2=results[[resname]]$predicted$M2), m2v_masks, data$M2, data$M3, data$M4, FALSE, FALSE)
        comoments_used <- c(comoments_used, "var")
      }
      loss <- as.numeric(loss)
      out[resname, paste0(loss_type,"_test_loss")] <- loss
      if (i > 0) {
        out[resname, paste0(loss_type,"_diff")] <- out[1, paste0(loss_type,"_test_loss")] - loss
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
      out[resname, paste0(loss_type,"_test_chisq")] <- data$meta_data$N * out[resname, paste0(loss_type,"_test_loss")]
      out[resname, paste0(loss_type,"_test_bic")] <- data$meta_data$N * out[resname, paste0(loss_type,"_test_loss")] + ncol(results[[resname]]$df) * log(data$meta_data$N)
    }
  }

  # Store N parameters loss type, use skew, use kurt
  for (resname in names(results)) {
    out[resname, "N_parameters"] <- ncol(results[[resname]]$df)
    for (i in c("loss_type", "use_skewness", "use_kurtosis")) {
      out[resname, i] <- results[[resname]]$info[[i]]
    }
  }
  if (!extensive_model_info) {
    for (i in c("loss_type", "use_skewness", "use_kurtosis")) {
      if (length(unique(out[, i]))==1) {out[,i] <- NULL}
    }
  } else {
    # Store optimizer, optim_iters and learning rates
    for (resname in names(results)) {
      for (n in seq_along(results[[resname]]$info$optimizers)) {
        out[resname, paste0("optim", n, "")] <- results[[resname]]$info$optimizers[n]
        out[resname, paste0("optim", n, "_iters")] <- results[[resname]]$info$optim_iters[n]
        out[resname, paste0("optim", n, "_lr")] <- results[[resname]]$info$learning_rate[n]
      }
    }
    # Store train_n version, use_bounds, optimizer, compute se, se type, device, low memory and runtimes
    for (resname in names(results)) {
      out[resname, 'train_n'] <- results[[resname]]$model$meta_data$n_obs
      for (i in c("version", "use_bounds", "compute_se", "se_type", "device", "low_memory")) {
        out[resname, i] <- results[[resname]]$info[[i]]
      }
      for (i in c("Preparation", "Optimizer", "SE", "Total")) {
        out[resname, paste0(i,"_runtime")] <- as.numeric(results[[resname]]$runtimes[[i]], units='secs')
      }
    }
  }
  return(out)
}