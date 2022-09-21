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
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)

    plot(mymodel)
}

