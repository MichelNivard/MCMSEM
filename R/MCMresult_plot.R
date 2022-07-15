plot.mcmresultclass <- function(res,layout=NULL) {
  # As of version 0.6.0 this is ported to plot.mcmmodelclass to prevent duplication of code
  return(plot(res$model,layout=NULL, use_values = TRUE))
}
