plot.mcmresultclass <- function(res,layout=NULL) {
  # As of version 0.6.0 this is ported to plot.mcmmodelclass to prevent duplication of code
  return(plot(res$model,layout=layout, use_values = TRUE))
}

plot.mcmresultsummaryclass <- function(summ,layout=NULL) {
  # Since users will inevitably try to plot the result summary, instead of the result, this will allow that
  return(plot(summ$result$model,layout=layout, use_values = TRUE))
}

