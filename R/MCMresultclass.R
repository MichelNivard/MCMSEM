mcmresultclass <- setRefClass("mcmresultclass",
                             fields=list(
                               df="data.frame",
                               loss="numeric",
                               history="list",
                               runtimes="list"
                             ))

mcmresultclass$methods(
  initialize=function(df, loss, history, runtimes) {
    .self$df <- df
    .self$loss <- loss
    .self$history <- history
    .self$runtimes <- runtimes
  },
  show=function(){
    cat("  MCM model Result\n")
    print(.self$df)
  },
  copy=function(){
    return(mcmresultclass(df, loss, history, runtimes))
  }
)

as.data.frame.mcmresultclass <- function(x) {
  return(x$df)
}