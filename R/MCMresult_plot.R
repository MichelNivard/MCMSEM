plot.mcmresultclass <- function(res) {
  # As of version 0.6.0 this is ported to plot.mcmmodelclass to prevent duplication of code
  return(plot(res$model, use_values = TRUE))
}
