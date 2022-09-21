\name{plot(mcmresultclass)}
\alias{plot(mcmresultclass)}
\title{Plot an MCM result}
\description{
  Plot qgraph-style path diagram representation of the MCM model results. This function is identical to plot(mcmmodelclass), but uses the fitted values to determine the width and color of the graph paths.
}
\usage{
plot(res, layout=NULL, ...)
}
\arguments{
  \item{res}{MCM result object}
  \item{layout}{qgraph style layout}
  \item{...}{Other arguments passed to qgraph}
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
    myresult <- MCMfit(model, data)

    plot(myresult)
}

