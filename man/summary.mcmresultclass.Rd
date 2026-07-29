\name{summary.mcmresultclass}
\alias{summary.mcmresultclass}
\title{Summary of MCMSEM results}
\description{
  Generate a kernel-aware summary of MCMSEM results. Dynamic summaries include
  transition paths, innovation cumulants, Gaussian residual covariance,
  spectral radius, stationarity, convergence, moment count, and nominal degrees
  of freedom.
}
\usage{
\method{summary}{mcmresultclass}(object, ...)
}
\arguments{
  \item{object}{MCM result object}
  \item{...}{Unused.}
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
