\name{as.data.frame.mcmresultsummaryclass}
\alias{as.data.frame.mcmresultsummaryclass}
\title{Extract dataframes from MCM result summary}
\description{
  Extract dataframes from MCM result summary
}
\usage{
\method{as.data.frame}{mcmresultsummaryclass}(x, row.names=NULL,
                                             optional=FALSE, ...,
                                             estimates='parameters')
}
\arguments{
  \item{x}{MCM result summary object}
  \item{row.names}{Unused; present for S3 compatibility.}
  \item{optional}{Unused; present for S3 compatibility.}
  \item{...}{Unused.}
  \item{estimates}{estimates to extract, should be one of "parameters", "variances", "skewness" or "kurtosis"}
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary, kernel="contemporaneous")
    my_result <- MCMfit(mymodel, mysummary, compute_se=FALSE,
                        optimizers="lbfgs", optim_iters=1, learning_rate=0.1)
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
