\name{as.data.frame(mcmresultsummaryclass)}
\alias{as.data.frame(mcmresultsummaryclass)}
\title{Extract dataframes from MCM result summary}
\description{
  Extract dataframes from MCM result summary
}
\usage{
as.data.frame(summ, estimates='parameters')
}
\arguments{
  \item{summ}{MCM result summary object}
  \item{estimates}{estimates to extract, should be one of "parameters", "variances", "skewness" or "kurtosis"}
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
