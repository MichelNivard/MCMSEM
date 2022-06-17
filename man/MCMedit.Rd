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
  \item{pointer}{Pointer to the location of the value(s) that should be changed. Either label of a matrix or one of c('bound', 'ubound', 'lbound', 'start').}
  \item{name}{Name (or location) of the parameter(s) that should be changed.}
  \item{value}{New value, this can be a string (to create a new parameter), or a number, to change specific values.}
}
\value{
  Returns a new MCMmodelclass object with changes applied
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)

    # free f1 variance, i.e. make it a parameter rather than hard value
    mymodel <- MCMedit(mymodel, "Fm", c(1, 2), "fm1")

    # Constrain b1_1 to zero
    mymodel <- MCMedit(mymodel, "A", "b1_1", 0)

    # Change lower bound of a1 to -.5
    mymodel <- MCMedit(mymodel, "lbound", "a1", -.5)

    # Change the upper bound of a1 to 2
    mymodel <- MCMedit(mymodel, "ubound", "a1", 2)

    # Set lower and upper bounds of a1 to c(-1, 1)
    mymodel <- MCMedit(mymodel, "bound", "a1", c(-1, 1))

    # Set lower and upper bounds of all b parameters to -1, 1
    mymodel <- MCMedit(mymodel, "bound", "b", c(-1, 1))

    # Set upper bound of all k parameters to 200
    mymodel <- MCMedit(mymodel, "ubound", "k", 200)

    # Set the starting value of b2 to 1.5
    mymodel <- MCMedit(mymodel, "start", "b2", 1.5)
    # Set the starting value of all a parameters to 1.0
    mymodel <- MCMedit(mymodel, "start", "a", 1.0)
}