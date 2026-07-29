\name{as.data.frame.mcmresultclass}
\alias{as.data.frame.mcmresultclass}
\title{Convert MCM result to dataframe}
\description{
  Convert MCM result to dataframe
}
\usage{
\method{as.data.frame}{mcmresultclass}(x, row.names=NULL,
                                      optional=FALSE, ...)
}
\arguments{
  \item{x}{MCM result object}
  \item{row.names}{Unused; present for S3 compatibility.}
  \item{optional}{Unused; present for S3 compatibility.}
  \item{...}{Unused.}
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary, kernel="contemporaneous")
    my_result <- MCMfit(mymodel, mysummary, compute_se=FALSE,
                        optimizers="lbfgs", optim_iters=1, learning_rate=0.1)
    my_result_df <- as.data.frame(my_result)
}
