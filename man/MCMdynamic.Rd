\name{MCMdynamic}
\alias{MCMmomentcount}
\alias{MCMdegreesoffreedom}
\alias{MCMdiagnostics}
\alias{MCMdynamicloss}
\title{Moment counts, identification diagnostics, and dynamic loss}
\description{
Inspect unique moments, independent-free-parameter counts, nominal degrees of
freedom, optional Jacobian diagnostics, dynamic stationarity, or dynamic loss.
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
  \item{jacobian}{Numerically evaluate the implied-moment Jacobian with respect
  to independent optimizer coordinates and report its rank and conditioning.}
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
number of independent free parameters. Fixed and derived parameters are
reported separately and do not consume degrees of freedom. Diagnostics include
the observed-moment count, optimizer-coordinate Jacobian, local rank, condition
number, near-null directions when detected, and the derived constraints used.
These local diagnostics do not prove global identification. Dynamic
identification also depends on useful
non-Gaussianity, independent innovations, fixed innovation scale, stationarity,
distinct source cumulant profiles, and no unmodelled non-Gaussian confounding.
The optional \code{residual_family="common_gamma"} makes one such confounder
explicit as a centered, variance-one gamma factor with signed loadings and a
positive shape. It contributes rank-one second-, third-, and fourth-order
cumulants directly to the stationary marginal distribution; it is not an
innovation propagated through \code{B}.
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
