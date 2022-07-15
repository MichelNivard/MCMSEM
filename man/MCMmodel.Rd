\name{MCMmodel}
\alias{MCMmodel}
\title{Create an MCMmodelclass object}
\description{
  Generate an MCMmodelclass object from the shape of observed data
}
\usage{
MCMmodel(data, n_confounding=1, constrained_a=TRUE, scale_data=TRUE, confounding_names=NULL)
}
\arguments{
  \item{data}{matrix or dataframe of the raw data the model should be built for}
  \item{n_confounding}{Number of latent confounding factors to be included.}
  \item{constrained_a}{Constrain a paths to be identical across observed variables (i.e. 1 parameter per confounder, as opposed to n_var parameters)}
  \item{scale_data}{Should the data be scaled prior to analysis?}
  \item{confounding_names}{Labels to be given to the latent confounding factors This is purely for the purposes of }
}
\value{
  Returns an MCMmodelclass object
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
}