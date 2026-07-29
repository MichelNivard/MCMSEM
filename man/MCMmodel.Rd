\name{MCMmodel}
\alias{MCMmodel}
\title{Create an MCMmodelclass object}
\description{
  Generate an MCMmodelclass object from the shape of observed data
}
\usage{
MCMmodel(data, n_latent=1, constrained_a=TRUE, scale_data=TRUE,
         weights=NULL, latent_names=NULL,
         causal_observed=TRUE, var_observed=TRUE, skew_observed=TRUE, kurt_observed=TRUE,
         causal_latent=FALSE, var_latent=FALSE, skew_latent=FALSE, kurt_latent=FALSE,
         kernel="contemporaneous", gaussian_residual=TRUE)
}
\arguments{
  \item{data}{Either an MCM data summary object, or a matrix or dataframe for which the model will be built}
  \item{n_latent}{Number of latent latent factors to be included.}
  \item{constrained_a}{Constrain a paths to be identical across observed variables (i.e. 1 parameter per confounder, as opposed to n_var parameters)}
  \item{scale_data}{Should the data be scaled prior to analysis? Only used if a matrix or dataframe is provided to data}
  \item{weights}{Weights to be used for weighted analysis. Should be vector with length equal to the number of rows in the raw data. Only used when a matrix or dataframe is provided to data}
  \item{latent_names}{Labels to be given to the latent factors This is only used for generation of table and figures.}
  \item{causal_observed}{Adds free parameters for causal paths between observed variables (default TRUE)}
  \item{var_observed}{Adds free parameters for variances of observed variables}
  \item{skew_observed}{Adds free parameters for skewness of observed variables}
  \item{kurt_observed}{Adds free parameters for kurtosis of observed variables}
  \item{causal_latent}{Adds free parameters for causal paths between latent factors}
  \item{var_latent}{Adds free parameters for variances of latent factors}
  \item{skew_latent}{Adds free parameters for skewness of latent factors}
  \item{kurt_latent}{Adds free parameters for kurtosis of latent factors}
  \item{kernel}{Exact kernel name: \code{"contemporaneous"} (the default) or
  \code{"dynamic"}. The deprecated alias \code{"static"} warns and maps to
  \code{"contemporaneous"}.}
  \item{gaussian_residual}{For a dynamic model, estimate a full Gaussian
  residual/random-intercept covariance using a Cholesky parameterization.}
}
\value{
  Returns an MCMmodelclass object
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary)

    # The dynamic kernel currently supports observed states only.
    dynamic_model <- MCMmodel(mysummary, kernel="dynamic")
}
