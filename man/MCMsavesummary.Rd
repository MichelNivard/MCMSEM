\name{MCMsavesummary}
\alias{MCMsavesummary}
\title{Save an MCM data summary}
\description{
  Save an MCM data summary disk
}
\usage{
MCMsavesummary(summaryobj, path, debug=FALSE)
}
\arguments{
  \item{summaryobj}{MCM data summary object}
  \item{path}{Path where to store the MCM dat summary (suffix '.mcmdata' is included automatically)}
  \item{debug}{Print more process information}
}

\examples{
    rawdata <- simulate_data(n=1000, asdataframe=TRUE)

    # Create summary
    mysummary <- MCMdatasummary(rawdata, prep_asymptotic_se=FALSE)

    # Save summary
    path <- tempfile()
    MCMsavesummary(mysummary, path)

    # Load summary
    my_other_summary <- MCMdatasummary(path=paste0(path, ".mcmdata"))
    unlink(paste0(path, ".mcmdata"))
}
