\name{simulate_data}
\alias{simulate_data}
\title{Simulate data for MCMSEM}
\description{
  Simulates data for a multi co-moment structural equation model
}
\usage{
simulate_data <- function(n=500000, a1=0.35, b1=0.3, b2=-.1, shape=4, df=10)
}
\arguments{
  \item{n}{Sample size}
  \item{a1}{Confouding path(s) (p -> x and p -> y)}
  \item{b1}{Causal path 1 (x -> y)}
  \item{b2}{Causal path 2 (y -> x)}
  \item{shape}{Skewness}
  \item{df}{Kurtosis}
}
\value{
  The function outputs a dataframe of shape [n, 2] matching the provided parameters. Note this is not exact, this simulated data will include sampling error.
}

\examples{
  simulate_data()
}
