\name{MCMmodel}
\alias{MCMmodel}
\title{Create an MCMmodelclass object}
\description{
  Generate an MCMmodelclass object from the shape of observed data
}
\usage{
MCMmodel(data, n_latent=1, constrained_a=TRUE, scale_data=TRUE, latent_names=NULL)
}
\arguments{
  \item{data}{matrix or dataframe of the raw data the model should be built for}
  \item{n_latent}{Number of latent latent factors to be included.}
  \item{constrained_a}{Constrain a paths to be identical across observed variables (i.e. 1 parameter per confounder, as opposed to n_var parameters)}
  \item{scale_data}{Should the data be scaled prior to analysis?}
  \item{latent_names}{Labels to be given to the latent factors This is only used for generation of table and figures.}
}
\value{
  Returns an MCMmodelclass object
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
}