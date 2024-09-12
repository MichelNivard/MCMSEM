\name{print(mcmmodelclass)}
\alias{print(mcmmodelclass)}
\title{Print MCM model matrices}
\description{
  Print MCM model matrices
}
\usage{
print(model, matrix=NULL)
}
\arguments{
  \item{model}{MCM model object}
  \item{matrix}{Label of the matrix to print, should be NULL or one of "A", "Fm", "S", "Sk" or "K". If NULL, A Fm and S are printed}
}

\examples{
    mysummary <- simulate_data()
    mymodel <- MCMmodel(mysummary)

    print(mymodel, matrix="A")
}
