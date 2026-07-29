\name{print.mcmmodelclass}
\alias{print.mcmmodelclass}
\title{Print MCM model matrices}
\description{
  Print MCM model matrices
}
\usage{
\method{print}{mcmmodelclass}(x, matrix=NULL, ...)
}
\arguments{
  \item{x}{MCM model object}
  \item{matrix}{Label of a matrix in the model. Contemporaneous matrices include
  "A", "Fm", "S", "Sk", and "K"; dynamic matrices include "B", "Tau",
  "Kappa", "L_G", and fixed "D2". If NULL, the principal matrices and kernel
  label are printed.}
  \item{...}{Unused.}
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary)

    print(mymodel, matrix="A")
}
