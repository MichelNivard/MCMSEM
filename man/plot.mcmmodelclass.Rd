\name{plot(mcmmodelclass)}
\alias{plot(mcmmodelclass)}
\title{Plot an MCM model}
\description{
  Plot qgraph-style path diagram representation of the MCM model
}
\usage{
plot(model, layout=NULL, ...)
}
\arguments{
  \item{model}{MCM model object}
  \item{layout}{qgraph style layout}
  \item{...}{Other arguments passed to qgraph}
}

\examples{
    mysummary <- simulate_data()
    mymodel <- MCMmodel(mysummary)

    plot(mymodel)
}

