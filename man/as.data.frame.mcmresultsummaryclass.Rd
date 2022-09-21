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
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
    my_result <- MCMfit(mymodel, data)
    summ <- summary(my_result)

    # Extract regression paramters and factor loadings:
    pars <- as.data.frame(summ, estimates='parameters')
    # Extract variances
    vars <- as.data.frame(summ, estimates='variances')
    # Extract skews
    skews <- as.data.frame(summ, estimates='skewness')
    # Extract kurts
    kurts <- as.data.frame(summ, estimates='kurtosis')
}
