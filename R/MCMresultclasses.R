## MCM Result class
mcmresultclass <- setRefClass("mcmresultclass",
                             fields=list(
                               df="data.frame",
                               model="mcmmodelclass",
                               loss="numeric",
                               history="list",
                               runtimes="list",
                               info="list"
                             ))

mcmresultclass$methods(
  initialize=function(df=data.frame(), model=mcmmodelclass(), loss=as.numeric(NA), history=list(), runtimes=list(), info=list()) {
    .self$df <- df
    .self$model <- model
    .self$loss <- loss
    .self$history <- history
    .self$runtimes <- runtimes
    .self$info <- info
  },
  show=function(){
    cat("  MCM model Result\n")
    print(.self$df)
  },
  copy=function(){
    return(mcmresultclass(.self$df, .self$model, .self$loss, .self$history, .self$runtimes, .self$info))
  }
)

as.data.frame.mcmresultclass <- function(x) {
  return(x$df)
}


## MCM Result Summary class
mcmresultsummaryclass <- setRefClass("mcmresultsummaryclass",
                             fields=list(
                               df="data.frame",
                               loss="numeric",
                               n_par="numeric",
                               n_obs="numeric",
                               chisq="numeric",
                               aic="numeric",
                               bic="numeric",
                               result="mcmresultclass"
                             ))


mcmresultsummaryclass$methods(
  initialize=function(df, loss, n_par, n_obs, chisq, aic, bic, result) {
    .self$df <- df
    .self$loss <- loss
    .self$n_par <- n_par
    .self$n_obs <- n_obs
    .self$chisq <- chisq
    .self$aic <- aic
    .self$bic <- bic
    .self$result <- result
  },
  show=function(){
    cat("|--------------------------------------|\n")
    cat(paste0("|  MCM Result Summary (MCMSEM v", .self$result$info$version, ")  |\n"))
    cat("|--------------------------------------|\n")
    cat(paste0("device         : ", .self$result$info$device, "\n"))
    cat(paste0("N phenotypes   : ", .self$result$model$meta_data$n_phenotypes, "\n"))
    cat(paste0("N latents      : ", .self$result$model$meta_data$n_confounding, "\n"))
    cat(paste0("N observations : ", .self$n_obs, "\n"))
    cat(paste0("N parameters   : ", .self$n_par, "\n"))
    if (.self$result$info$compute_se) {cat(paste0("SE type        : ", .self$result$info$se_type, "\n"))}
    cat("\n")
    cat("Fit statistics\n")
    cat(paste0("loss  : ", .self$loss, "\n"))
    cat(paste0("chisq : ", .self$chisq, "\n"))
    cat(paste0("aic   : ", .self$aic, "\n"))
    cat(paste0("bic   : ", .self$bic, "\n"))
    cat("\n")
    cat("Parameter summary\n")
    if (nrow(.self$df) > 16) {
      print(.self$df[1:15, ])
      cat("...\n")
    } else {
      print(.self$df)
    }

  },
  copy=function(){
    return(mcmresultclass(.self$df, .self$loss, .self$n_par, .self$n_obs, .self$chisq, .self$aic, .self$bic, .self$result$copy()))
  }
)

as.data.frame.mcmresultsummaryclass <- function(x) {
  return(x$df)
}