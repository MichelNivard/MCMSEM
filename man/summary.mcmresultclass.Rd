\name{summary(mcmresultclass)}
\alias{summary(mcmresultclass)}
\title{Summary of MCMSEM results}
\description{
  Generate a summary of MCMSEM results
}
\usage{
summary(res)
}
\arguments{
  \item{res}{MCM result object}
}

\examples{
    mysummary <- simulate_data()
    mymodel <- MCMmodel(mysummary)
    my_result <- MCMfit(mymodel, mysummary)
    resultsummary <- summary(my_result)
    # Extract regression paramters and factor loadings:
    pars <- as.data.frame(resultsummary, estimates='parameters')
    # Extract variances
    vars <- as.data.frame(resultsummary, estimates='variances')
    # Extract skews
    skews <- as.data.frame(resultsummary, estimates='skewness')
    # Extract kurts
    kurts <- as.data.frame(resultsummary, estimates='kurtosis')
}
