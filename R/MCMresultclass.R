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
  initialize=function(df, model, loss, history, runtimes, info) {
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