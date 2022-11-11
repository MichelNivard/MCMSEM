\name{MCMdatasummary}
\alias{MCMdatasummary}
\title{Generate an MCM data summary}
\description{
  Create a new MCM data summary, or load an existing one from disk
}
\usage{
MCMdatasummary(data=NULL, path=NULL, weights=NULL, scale_data=TRUE, prep_asymptotic_se=TRUE, use_skewness=TRUE,
               use_kurtosis=TRUE, debug=FALSE, low_memory=FALSE)
}
\arguments{
  \item{data}{Raw data from which to create a new MCM data summary object (either data or path must be provided)}
  \item{path}{Path to an existing data summary object (either data or path must be provided)}
  \item{weights}{Weights to be used for generating co-moment matrices for weighted analysis}
  \item{scale_data}{Should the data be scaled prior to analysis?}
  \item{prep_asymptotic_se}{Prepare the required matrix for calculation of asymptotic standard errors}
  \item{use_skewness}{Include matrices for a model using skewness}
  \item{use_kurtosis}{Include matrices for a model using kurtosis}
  \item{debug}{Print detailed process updates}
  \item{low_memory}{Integer to be set when running into memory issues ranging from 0 (fast, no memory saving) to 4 (slow, highest memory saving)}
}
\value{
  An MCM data summary object which can be sent to collaborators instead of the raw data. This object does not include individual records, but is everything that is required to run all other functions in MCMSEM (except for bootstrapped standard errors)
}

\examples{
    rawdata <- simulate_data(asdataframe=TRUE)

    # Create summary
    mysummary <- MCMdatasummary(rawdata)

    # Save summary
    MCMsavesummary(mysummary, "path/to/summary")

    # Load summary
    my_other_summary <- MCMdatasummary(path="path/to/summary.mcmdata")
}
