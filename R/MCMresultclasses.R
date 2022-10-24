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
                               predicted="list"
                             ))

mcmresultclass$methods(
  initialize=function(df=data.frame(), model=mcmmodelclass(), loss=as.numeric(NA), gradients=mcmmultigradienthistoryclass(), history=list(), runtimes=list(), info=list(),
                      observed=list(), predicted=list()) {
    .self$df <- df
    .self$model <- model
    .self$loss <- loss
    .self$gradients <- gradients
    .self$history <- history
    .self$runtimes <- runtimes
    .self$info <- info
    .self$observed <- observed
    .self$predicted <- predicted
  },
  show=function(){
    cat("  MCM model Result\n")
    print(.self$df)
  },
  copy=function(){
    return(mcmresultclass(.self$df, .self$model$copy(), .self$loss, .self$gradients, .self$history, .self$runtimes, .self$info, .self$observed, .self$predicted))
  }
)

as.data.frame.mcmresultclass <- function(x) {
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
    cat(paste0("device         : ", .self$result$info$device, "\n"))
    cat(paste0("N phenotypes   : ", .self$result$model$meta_data$n_phenotypes, "\n"))
    cat(paste0("N latents      : ", .self$result$model$meta_data$n_latent, "\n"))
    cat(paste0("N observations : ", .self$n_obs, "\n"))
    cat(paste0("N parameters   : ", .self$n_par, "\n"))
    if (.self$result$info$compute_se) {cat(paste0("SE type        : ", .self$result$info$se_type, "\n"))}
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
    return(mcmresultclass(.self$parameters, .self$variances, .self$skewness, .self$kurtosis, .self$loss, .self$n_par, .self$n_obs, .self$chisq, .self$bic, .self$result$copy()))
  }
)

as.data.frame.mcmresultsummaryclass <- function(summ, estimates="parameters") {
  if (estimates %in% c("parameters", "variances", "skewness", "kurtosis")) {
    return(summ[[estimates]])
  } else {
    stop('estimates argument should be one of ("parameters", "variances", "skewness", "kurtosis")')
  }
}