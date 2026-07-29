## MCM Result class
mcmresultclass <- setRefClass("mcmresultclass",
                             fields=list(
                               df="data.frame",
                               model="mcmmodelclass",
                               loss="numeric",
                               gradients="mcmmultigradienthistoryclass",
                               history="list",
                               runtimes="list",
                               info="list",
                               observed="list",
                               predicted="list",
                               kernel="character",
                               B="matrix",
                               transition_parameters="data.frame",
                               innovation_variances="vector",
                               innovation_third="vector",
                               innovation_fourth="vector",
                               Psi_G="matrix",
                               spectral_radius="numeric",
                               stationary="logical",
                               convergence="list",
                               degrees_of_freedom="numeric",
                               n_moments="numeric",
                               start_diagnostics="data.frame",
                               residuals="list",
                               dynamic="list"
                             ))

mcmresultclass$methods(
  initialize=function(df=data.frame(), model=mcmmodelclass(), loss=as.numeric(NA), gradients=mcmmultigradienthistoryclass(), history=list(), runtimes=list(), info=list(),
                      observed=list(), predicted=list(), kernel="contemporaneous",
                      B=matrix(), transition_parameters=data.frame(),
                      innovation_variances=numeric(), innovation_third=numeric(),
                      innovation_fourth=numeric(), Psi_G=matrix(),
                      spectral_radius=as.numeric(NA), stationary=FALSE,
                      convergence=list(), degrees_of_freedom=as.numeric(NA),
                      n_moments=as.numeric(NA), start_diagnostics=data.frame(),
                      residuals=list(), dynamic=list()) {
    .self$df <- df
    .self$model <- model
    .self$loss <- loss
    .self$gradients <- gradients
    .self$history <- history
    .self$runtimes <- runtimes
    .self$info <- info
    .self$observed <- observed
    .self$predicted <- predicted
    .self$kernel <- kernel
    .self$B <- B
    .self$transition_parameters <- transition_parameters
    .self$innovation_variances <- innovation_variances
    .self$innovation_third <- innovation_third
    .self$innovation_fourth <- innovation_fourth
    .self$Psi_G <- Psi_G
    .self$spectral_radius <- spectral_radius
    .self$stationary <- stationary
    .self$convergence <- convergence
    .self$degrees_of_freedom <- degrees_of_freedom
    .self$n_moments <- n_moments
    .self$start_diagnostics <- start_diagnostics
    .self$residuals <- residuals
    .self$dynamic <- dynamic
  },
  show=function(){
    cat("  MCM model Result\n")
    cat("Kernel: ", .kernel_label(.result_kernel(.self)), "\n", sep = "")
    print(.self$df)
  },
  copy=function(){
    field_or <- function(name, default) {
      tryCatch(.self[[name]], error = function(e) default)
    }
    return(mcmresultclass(
      .self$df, .self$model$copy(), .self$loss, .self$gradients,
      .self$history, .self$runtimes, .self$info, .self$observed,
      .self$predicted, .result_kernel(.self), field_or("B", matrix()),
      field_or("transition_parameters", data.frame()),
      field_or("innovation_variances", numeric()),
      field_or("innovation_third", numeric()),
      field_or("innovation_fourth", numeric()),
      field_or("Psi_G", matrix()),
      field_or("spectral_radius", as.numeric(NA)),
      field_or("stationary", FALSE), field_or("convergence", list()),
      field_or("degrees_of_freedom", as.numeric(NA)),
      field_or("n_moments", as.numeric(NA)),
      field_or("start_diagnostics", data.frame()),
      field_or("residuals", list()), field_or("dynamic", list())
    ))
  }
)

as.data.frame.mcmresultclass <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  return(x$df)
}


## MCM Result Summary class
mcmresultsummaryclass <- setRefClass("mcmresultsummaryclass",
                             fields=list(
                               parameters="data.frame",
                               variances="data.frame",
                               skewness="data.frame",
                               kurtosis="data.frame",
                               loss="numeric",
                               n_par="numeric",
                               n_obs="numeric",
                               chisq="numeric",
                               bic="numeric",
                               result="mcmresultclass"
                             ))


mcmresultsummaryclass$methods(
  initialize=function(parameters, variances, skewness, kurtosis, loss, n_par, n_obs, chisq, bic, result) {
    .self$parameters <- parameters
    .self$variances <- variances
    .self$skewness <- skewness
    .self$kurtosis <- kurtosis
    .self$loss <- loss
    .self$n_par <- n_par
    .self$n_obs <- n_obs
    .self$chisq <- chisq
    .self$bic <- bic
    .self$result <- result
  },
  show=function(){
    cat("|--------------------------------------|\n")
    cat(paste0("| MCM Result Summary (MCMSEM v", .self$result$info$version, ")",paste0(rep(" ",8-nchar(.self$result$info$version)), collapse=''),"|\n"))
    cat("|--------------------------------------|\n")
    cat("Kernel         : ", .kernel_label(.result_kernel(.self$result)), "\n", sep = "")
    cat(paste0("device         : ", .self$result$info$device, "\n"))
    cat(paste0("N phenotypes   : ", .self$result$model$meta_data$n_phenotypes, "\n"))
    cat(paste0("N latents      : ", .self$result$model$meta_data$n_latent, "\n"))
    cat(paste0("N observations : ", .self$n_obs, "\n"))
    cat(paste0("N parameters   : ", .self$n_par, "\n"))
    if (identical(.result_kernel(.self$result), "dynamic")) {
      cat(paste0("N moments      : ", .self$result$n_moments, "\n"))
      cat(paste0("Nominal df     : ", .self$result$degrees_of_freedom, "\n"))
      cat(paste0("Spectral radius: ", signif(.self$result$spectral_radius, 6), "\n"))
      cat(paste0("Stationary     : ", .self$result$stationary, "\n"))
      cat(paste0("Best start     : ", .self$result$convergence$best_start,
                 "/", .self$result$info$n_starts, "\n"))
      cat(paste0("Convergence    : ", .self$result$convergence$message, "\n"))
      cat(paste0("Moment weights : ", .self$result$info$moment_weighting, "\n"))
    }
    if (.self$result$info$compute_se) {
      cat(paste0("SE type        : ", .self$result$info$se_type,
                 " (", .self$result$info$se_correction, ")\n"))
      cat(paste0("Jacobian rank  : ", .self$result$info$jacobian_rank,
                 "/", length(.self$result$model$param_values), "\n"))
      cat(paste0("Info condition : ",
                 signif(.self$result$info$information_condition, 6), "\n"))
    }
    cat("\n")
    cat("Fit statistics\n")
    cat(paste0("loss  : ", .self$loss, "\n"))
    cat(paste0("chisq : ", .self$chisq, "\n"))
    cat(paste0("BIC   : ", .self$bic, "\n"))
    for (summ in c("parameters", "variances", "skewness", "kurtosis")) {
      if (nrow(.self[[summ]]) > 0) {
        cat("\n")
        cat(paste0(toupper(substr(summ, 1, 1)), substr(summ, 2, nchar(summ))), "summary\n")
        if (nrow(.self[[summ]]) > 16) {
          print(.self[[summ]][1:15, ])
          cat(paste0("... Print capped at 15 rows, use: as.data.frame([summary_object], estimates='",summ,"')\n"))
        } else {
          print(.self[[summ]])
        }
      }
    }
  },
  copy=function(){
    return(mcmresultsummaryclass(
      .self$parameters, .self$variances, .self$skewness,
      .self$kurtosis, .self$loss, .self$n_par, .self$n_obs,
      .self$chisq, .self$bic, .self$result$copy()
    ))
  }
)

as.data.frame.mcmresultsummaryclass <- function(x, row.names = NULL,
                                                optional = FALSE, ...,
                                                estimates="parameters") {
  if (estimates %in% c("parameters", "variances", "skewness", "kurtosis")) {
    return(x[[estimates]])
  } else {
    stop('estimates argument should be one of ("parameters", "variances", "skewness", "kurtosis")')
  }
}
