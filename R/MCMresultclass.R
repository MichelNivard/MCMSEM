mcmresultclass <- setRefClass("mcmresultclass",
                             fields=list(
                               df="data.frame",
                               loss="numeric",
                               history="list"
                             ))
mcmresultclass$methods(
  initialize=function(df, loss, history) {
    .self$df <- df
    .self$loss <- loss
    .self$history <- history
  },
  show=function(){
    cat("  MCM model Result\n")
    print(.self$df)
  }
)

as.data.frame.mcmresultclass <- function(x) {
  return(x$df)
}