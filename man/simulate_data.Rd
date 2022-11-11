\name{simulate_data}
\alias{simulate_data}
\title{Simulate data for MCMSEM}
\description{
  Simulates data for a multi co-moment structural equation model
}
\usage{
simulate_data <- function(n=500000, a=matrix(c(0.35, 0.35), ncol=1), b=matrix(c(1,-.1, .3,1), 2, 2),
                          shape=c(4, 0), df=c(0, 10),
                          asdataframe=FALSE, ...)
}
\arguments{
  \item{n}{Sample size}
  \item{a}{Matrix of confouding path(s). Should be of shape (n_variables, n_latents)}
  \item{b}{Matrix of causal paths. Paths should be on the off-diagonal, with a diagonal of 1. Should be of shape (n_variables, n_variables)}
  \item{shape}{Vector of shapes used for generating skewness. 0 indicates no skewness. Should be of length n_variables}
  \item{df}{Vector of dfs used for kurtosis. 0 indicates no kurtosis. Should be of length n_variables}
  \item{asdataframe}{Return the raw simulated data as data frame. If FALSE simulate_data() will return an MCMdatasummary class object}
  \item{...}{Additional arguments passed to MCMdatasummary() when asdataframe is FALSE}
}
\value{
  The function outputs a dataframe of shape (n, n_variables) matching the provided parameters. Note this is not exact, this simulated data will include sampling error.
}

\examples{
  simulate_data()
}
