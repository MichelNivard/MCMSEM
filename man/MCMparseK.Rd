\name{MCMparseK}
\alias{MCMparseK}
\title{Parse K matrix}
\description{
  Calculate the full product of S parameters in the K matrix of the given model.
}
\usage{
MCMparseK(model, sort=TRUE, shorten=TRUE, add_one=TRUE, print=FALSE)
}
\arguments{
  \item{model}{MCMmodelclass object, note that models with more than 5 observed variables will likely take a while to parse}
  \item{sort}{Sort parameters per multiplication operation, e.g. "s2+s2*s1+s1" -> "s2+s1*s2+s1"}
  \item{shorten}{Shorten the resulting formula where possible, e.g. "s2+s2+s2+s1*s1" -> "3*s2+s1^2", note that "sqrt(s1)*sqrt(s1)" will always be shortened to "s1" regardless of this setting}
  \item{add_one}{Add the "1*" prefix to cells with just a single parameter, e.g. "s1" -> "1*s1", to reproduce the equation in the original publication}
  \item{print}{Print the parsed K matrix instead of returning it, this will (attempt to) format K such that there's an equal number of columns printed on each row, and print without quotes to make the table as readable as possible}
}
\value{
  Returns a representation of the K parameters, with the product of the parameters in the S matrix and the free parameters in the K matrix.
}

\examples{
    mysummary <- simulate_data()
    mymodel <- MCMmodel(mysummary)

    print(MCMparseK(mymodel))
}
