mcmresultclass <- setRefClass("mcmresultclass",
                             fields=list(
                               df="data.frame",
                               model="mcmmodelclass",
                               loss="numeric",
                               history="list",
                               runtimes="list"
                             ))

mcmresultclass$methods(
  initialize=function(df, model, loss, history, runtimes) {
    .self$df <- df
    .self$model <- model
    .self$loss <- loss
    .self$history <- history
    .self$runtimes <- runtimes
  },
  show=function(){
    cat("  MCM model Result\n")
    print(.self$df)
  },
  copy=function(){
    return(mcmresultclass(.self$df, .self$model, .self$loss, .self$history, .self$runtimes))
  }
)

as.data.frame.mcmresultclass <- function(x) {
  return(x$df)
}