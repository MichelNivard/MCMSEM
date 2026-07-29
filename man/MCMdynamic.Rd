\name{MCMdynamic}
\alias{MCMmomentcount}
\alias{MCMdegreesoffreedom}
\alias{MCMdiagnostics}
\alias{MCMdynamicloss}
\title{Dynamic-kernel counts, diagnostics, and loss}
\description{
Inspect the unique moments, nominal degrees of freedom, stationarity, optional
Jacobian rank, or loss at supplied parameters for Stationary Dynamic MCMSEM.
}
\usage{
MCMmomentcount(x, use_skewness=TRUE, use_kurtosis=TRUE)

MCMdegreesoffreedom(object, use_skewness=TRUE, use_kurtosis=TRUE)

MCMdiagnostics(object, jacobian=FALSE,
               tolerance=sqrt(.Machine$double.eps))

MCMdynamicloss(model, data, parameters=NULL, loss_type="mse",
               moment_weighting="identity", weight_ridge=1e-8)
}
\arguments{
  \item{x}{A positive variable count, MCM data summary, model, result, matrix,
  or data frame.}
  \item{use_skewness}{Include unique third-order moments.}
  \item{use_kurtosis}{Include unique fourth-order moments.}
  \item{object}{An MCM model or result.}
  \item{jacobian}{For a dynamic model, numerically evaluate the implied-moment
  Jacobian and report its rank.}
  \item{tolerance}{Rank tolerance passed to \code{qr()}.}
  \item{model}{A model made with \code{kernel="dynamic"}.}
  \item{data}{An MCM data summary or raw numeric data.}
  \item{parameters}{Optional parameter vector in \code{model$param_names}
  order. The model's current values are used by default.}
  \item{loss_type}{One of \code{"mse"}, \code{"l1"}, or
  \code{"smooth_l1"}.}
  \item{moment_weighting}{One of \code{"identity"}, \code{"diagonal"}, or
  \code{"full"}. Non-identity weights require a prepared moment covariance and
  \code{loss_type="mse"}.}
  \item{weight_ridge}{Non-negative relative floor used when constructing
  diagonal or full moment weights.}
}
\value{
\code{MCMmomentcount()} returns named counts. \code{MCMdegreesoffreedom()} and
\code{MCMdiagnostics()} return diagnostic lists. \code{MCMdynamicloss()} returns
a numeric loss.
}
\details{
Nominal degrees of freedom are the number of unique observed moments minus the
number of free parameters. This count, and a local numerical Jacobian rank, do
not prove global identification. Dynamic identification also depends on useful
non-Gaussianity, independent innovations, fixed innovation scale, stationarity,
distinct source cumulant profiles, and no unmodelled non-Gaussian confounding.
}
\references{
Tamimy, Zenab, Elsje van Bergen, Matthijs D. van der Zee, Conor V. Dolan, and
Michel G. Nivard (2022). \emph{Multi Co-Moment Structural Equation Models:
Discovering Direction of Causality in the Presence of Confounding}. SocArXiv.
\doi{10.31235/osf.io/ynam2}.

Recke, Cecilie Olesen, Sarah Lumpp, Nataliia Kushnerchuk, Janike Oldekop,
Jiayi Li, Jane Ivy Coons, and Elina Robeva (2026). \emph{Identifiability in
Graphical Discrete Lyapunov Models}. arXiv preprint arXiv:2601.21818.
\url{https://arxiv.org/abs/2601.21818}.
}
\examples{
  ds <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
  model <- MCMmodel(ds, kernel="dynamic")
  MCMmomentcount(model)
  MCMdegreesoffreedom(model)
}
