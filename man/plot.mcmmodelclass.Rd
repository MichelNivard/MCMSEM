\name{plot.mcmmodelclass}
\alias{plot.mcmmodelclass}
\title{Plot an MCM model}
\description{
  Plot qgraph-style path diagram representation of the MCM model
}
\usage{
\method{plot}{mcmmodelclass}(x, layout=NULL, use_values=FALSE, ...)
}
\arguments{
  \item{x}{MCM model object}
  \item{layout}{qgraph style layout}
  \item{use_values}{Use fitted numeric path values rather than labels.}
  \item{...}{Other arguments passed to qgraph}
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary, kernel="contemporaneous")

    plot(mymodel)
}
