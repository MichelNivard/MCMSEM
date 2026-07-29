\name{plot.mcmresultclass}
\alias{plot.mcmresultclass}
\title{Plot an MCM result}
\description{
  Plot qgraph-style path diagram representation of the MCM model results. This function is identical to plot(mcmmodelclass), but uses the fitted values to determine the width and color of the graph paths.
}
\usage{
\method{plot}{mcmresultclass}(x, layout=NULL, ...)
}
\arguments{
  \item{x}{MCM result object}
  \item{layout}{qgraph style layout}
  \item{...}{Other arguments passed to qgraph}
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary)
    myresult <- MCMfit(mymodel, mysummary, compute_se=FALSE,
                       optimizers="lbfgs", optim_iters=1, learning_rate=0.1)

    plot(myresult)
}
