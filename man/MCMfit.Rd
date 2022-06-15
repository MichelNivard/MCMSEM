\name{MCMfit}
\alias{MCMfit}
\title{Fit an MCM model}
\description{
  Fits an MCM sturctural equation model as defined in the MCMmodelclass object
}
\usage{
MCMfit(model, data, compute_se=TRUE, se_type='asymptotic', bootstrap_iter=200,bootstrap_chunks=1000)
}
\arguments{
  \item{model}{MCMmodelclass object}
  \item{data}{The matrix or dataframe that was used to generate the MCMmodelclass object}
  \item{compute_se}{Compute SE for the parameters}
  \item{se_type}{How the SE should be calculated. Should be one of 'asymptotic': ..., 'two-step':... or 'one-step':}
  \item{bootstrap_iter}{Number of bootstrap iterations}
  \item{bootstrap_chunk}{Chunksize of two-step bootstrap. Lower numbers will lead to better performance but may lead to less accurate estimates.}
}
\value{
  Dataframe with parameter estimates, and SEs if compute_se is set to TRUE
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
    my_result <- MCMfit(mymodel)
}