\name{MCMSEM-classes}
\alias{mcmdataclass}
\alias{mcmdataclass-class}
\alias{mcmmodelclass}
\alias{mcmmodelclass-class}
\alias{mcmstartvaluesclass}
\alias{mcmstartvaluesclass-class}
\alias{mcmresultclass}
\alias{mcmresultclass-class}
\alias{mcmresultsummaryclass}
\alias{mcmresultsummaryclass-class}
\alias{mcmgradienthistoryclass}
\alias{mcmgradienthistoryclass-class}
\alias{mcmmultigradienthistoryclass}
\alias{mcmmultigradienthistoryclass-class}
\alias{MCMSEMversion}
\alias{summary.mcmmodelclass}
\alias{plot.mcmresultsummaryclass}
\alias{abs.mcmgradienthistoryclass}
\alias{log.mcmgradienthistoryclass}
\alias{log10.mcmgradienthistoryclass}
\alias{sqrt.mcmgradienthistoryclass}
\alias{min.mcmgradienthistoryclass}
\alias{max.mcmgradienthistoryclass}
\alias{summary.mcmgradienthistoryclass}
\alias{plot.mcmgradienthistoryclass}
\alias{as.data.frame.mcmgradienthistoryclass}
\alias{as.matrix.mcmgradienthistoryclass}
\alias{abs.mcmmultigradienthistoryclass}
\alias{log.mcmmultigradienthistoryclass}
\alias{log10.mcmmultigradienthistoryclass}
\alias{sqrt.mcmmultigradienthistoryclass}
\alias{min.mcmmultigradienthistoryclass}
\alias{max.mcmmultigradienthistoryclass}
\alias{summary.mcmmultigradienthistoryclass}
\alias{plot.mcmmultigradienthistoryclass}
\alias{as.data.frame.mcmmultigradienthistoryclass}
\alias{as.matrix.mcmmultigradienthistoryclass}
\title{MCMSEM support classes and methods}
\description{
Reference classes used for MCM data summaries, models, fitted results, result
summaries, start values, and gradient histories, together with their supporting
S3 methods. These objects are normally created by \code{MCMdatasummary()},
\code{MCMmodel()}, \code{MCMfit()}, and \code{summary()} rather than by calling
the class generators directly.
}
\details{
\code{MCMSEMversion} contains the package's internal result-version label.
Gradient-history transformation and plotting methods operate on histories
created when \code{MCMfit(monitor_grads=TRUE)} is used.
}
