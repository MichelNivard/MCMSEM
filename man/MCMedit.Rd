\name{MCMedit}
\alias{MCMedit}
\title{Edit existing MCM models}
\description{
  Edit parameters, starting values, or bounds of existing models
}
\usage{
MCMedit(model, pointer, name, value)
}
\arguments{
  \item{model}{Existing MCMmodelclass object}
  \item{pointer}{Pointer to the location of the value(s) that should be changed.
  Either a matrix label or one of c('bound', 'ubound', 'lbound', 'start').
  Dynamic matrices are \code{B}, \code{Tau}, \code{Kappa}, and \code{L_G}.}
  \item{name}{Name (or location) of the parameter(s) that should be changed.}
  \item{value}{New value, this can be a string (to create a new parameter), or a number, to change specific values.}
}
\value{
  Returns a new MCMmodelclass object with changes applied
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary, kernel="dynamic")

    # Rename a cross-lagged transition and edit its start and bounds.
    mymodel <- MCMedit(mymodel, "B", c(1, 2), "Y_lag_to_X")
    mymodel <- MCMedit(mymodel, "start", "Y_lag_to_X", 0.1)
    mymodel <- MCMedit(mymodel, "bound", "Y_lag_to_X", c(-0.8, 0.8))
}
