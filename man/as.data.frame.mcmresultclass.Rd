\name{as.data.frame(mcmresultclass)}
\alias{as.data.frame(mcmresultclass)}
\title{Convert MCM result to dataframe}
\description{
  Convert MCM result to dataframe
}
\usage{
as.data.frame(res)
}
\arguments{
  \item{res}{MCM result object}
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
    my_result <- MCMfit(mymodel, data)
    my_result_df <- as.data.frame(my_result)
}
