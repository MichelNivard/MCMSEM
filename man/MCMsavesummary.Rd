\name{MCMsavesummary}
\alias{MCMsavesummary}
\title{Save an MCM data summary}
\description{
  Save an MCM data summary disk
}
\usage{
MCMsavesummary(summaryobj, path)
}
\arguments{
  \item{summaryobj}{MCM data summary object}
  \item{path}{Path where to store the MCM dat summary (suffix '.mcmdata' is included automatically)}
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
