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
