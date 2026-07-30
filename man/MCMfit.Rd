\name{MCMfit}
\alias{MCMfit}
\title{Fit an MCM model}
\description{
  Fits an MCM sturctural equation model as defined in the MCMmodelclass object
}
\usage{
MCMfit(mcmmodel, data, weights=NULL, compute_se=TRUE,
       se_type='asymptotic', optimizers=c("rprop", "lbfgs"),
       optim_iters=c(50, 12), loss_type='mse', bootstrap_iter=200,
       bootstrap_chunks=1000, learning_rate=c(0.02, 1),
       use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE,
       device=NULL, device_se=NULL, low_memory=FALSE,
       outofbounds_penalty=1, monitor_grads=FALSE, jacobian_method='simple', debug=FALSE,
       n_starts=1L, seed=NULL, stationarity_limit=0.995,
       stationarity_penalty=1e6, verbose=debug,
       moment_weighting="identity", se_correction="auto",
       weight_ridge=1e-8)
}
\arguments{
  \item{mcmmodel}{MCMmodelclass object}
  \item{data}{Either an MCM data summary object, or a matrix or dataframe that was used to generate the MCMmodelclass object}
  \item{weights}{Weights to be used for weighted analysis. Should be vector with length equal to the number of rows in the raw data. Only used when a matrix or dataframe is provided to data}
  \item{compute_se}{Compute SE for the parameters}
  \item{se_type}{How the SE should be calculated. Should be one of 'asymptotic': asymptotic approximation, 'two-step': chunked bootstrap or 'one-step': full standard bootstrap}
  \item{optimizers}{Vector of optimizers to be used sequentially. The following optimizers are supported: c('rprop', 'sgd' ,'rmsprop', 'asgd', 'adam', 'adagrad', 'adadelta', 'lbfgs').}
  \item{optim_iters}{Number of iterations for each optimizer. Should be of length 1 or length(optimizers)}
  \item{learning_rate}{Learning rates for the optimizers defined in optimizers. Should be of length 1 or length(optimizers)}
  \item{loss_type}{How to compute loss, should be one of 'mse' (mean squared errors), or 'smooth_l1' (smooth L1).}
  \item{bootstrap_iter}{Number of bootstrap iterations}
  \item{bootstrap_chunks}{Number of chunks in the two-step bootstrap. Lower numbers can improve performance but may reduce accuracy.}
  \item{use_bounds}{Use parameter bounds defined in the model, set this to FALSE in larger models to improve compute times}
  \item{use_skewness}{Use skewness and co-skewness to estimate parameters}
  \item{use_kurtosis}{Use kurtosis and co-kurtosis to estimate parameters}
  \item{device}{torch device object, can be used to run optimization on GPU; device=torch_device('cuda'), see "MCMSEM on GPU" in the README, defaults to using torch_device('cpu').}
  \item{device_se}{torch device object used during calculation of SE (defaults to the same as device). Primarily used for changing device_se to a CPU device.}
  \item{low_memory}{Forces a more memory-friendly version of the optimization at the cost of runtime (note this cost gets VERY significant with large (>30 variables) models). Only set this option to TRUE when running into CUDA out of Memory errors.}
  \item{outofbounds_penalty}{loss penalty scaling when parameters are out of bounds: loss *  2 ^ (distance_from_bounds * outofbounds_penalty)}
  \item{monitor_grads}{Monitor the gradients during optimization and return gradients at each optimization. This will aditionally cause non-LBFGS optimization to stop early when NaN gradients are encountered, and returns the last non-NaN result instead.}
  \item{jacobian_method}{Method used for the `numDeriv::jacobian()` function in asymptotic SE calculation}
  \item{debug}{Print more process information, including loss at every iteration}
  \item{n_starts}{For dynamic models, total number of screened optimization
  starts. Failed starts do not terminate the fit.}
  \item{seed}{Optional seed for reproducible dynamic starts and torch
  optimization.}
  \item{stationarity_limit}{Dynamic-kernel spectral-radius limit, strictly
  between zero and one.}
  \item{stationarity_penalty}{Positive smooth penalty multiplier for dynamic
  parameter values beyond \code{stationarity_limit}.}
  \item{verbose}{For dynamic fits, print start number, initial/final loss,
  convergence, spectral radius, and residual covariance. Defaults to
  \code{debug}.}
  \item{moment_weighting}{Dynamic-kernel WLS weight: \code{"identity"},
  \code{"diagonal"}, or \code{"full"}. Diagonal and full weighting require a
  summary prepared for asymptotic SEs.}
  \item{se_correction}{Dynamic asymptotic covariance: \code{"auto"} selects
  efficient inverse-information SEs for full WLS and robust sandwich SEs
  otherwise; \code{"robust"} always uses the sandwich; \code{"model_based"}
  requires full WLS.}
  \item{weight_ridge}{Non-negative relative eigenvalue or diagonal floor used
  when forming diagonal or full WLS moment weights. The parameter information
  matrix itself is inverted without this statistical ridge.}
}
\value{
  An MCM result object. For constrained models, \code{df} and
  \code{parameter_table} report free, fixed, and derived estimates;
  \code{parameter_vcov} retains their joint reporting-scale covariance and
  \code{free_parameter_vcov} its independent-free-parameter block. The object
  can be converted to a data frame with \code{as.data.frame(result)}.
}
\details{
Dynamic models require skewness and kurtosis. Their asymptotic calculations use
the Jacobian of the implied raw central moments and the mean-corrected sampling
covariance prepared by \code{MCMdatasummary()}. Identity and diagonal weights
use the robust sandwich covariance by default; full WLS uses the efficient
inverse-information covariance. Residual covariance SEs are delta-method
transformations of the underlying Gaussian Cholesky or common-gamma loading
parameters. Analysis weights are not yet
supported for dynamic WLS or SEs. Dynamic results retain all-start diagnostics
and the best finite, in-bounds stationary solution.

For \code{residual_family="common_gamma"}, the fitted factor loadings and shape
are in \code{dynamic$common_gamma}; \code{dynamic$residual_M2},
\code{residual_M3}, and \code{residual_K4} retain its cumulant contribution.
Shape-implied skewness and excess-kurtosis SEs also use the delta method.
The legacy top-level \code{Psi_G} field remains reserved for Gaussian residuals.

Only independent free parameters enter either optimizer. Positive and bounded
free parameters are optimized on unconstrained internal coordinates and
reported on their natural scale. Estimator covariance is first computed in
those coordinates, then transformed to natural free and derived parameters by
the delta method. Fixed-parameter SEs are zero when covariance is available.
For dynamic results, \code{dynamic$asymptotic$vcov_robust} and
\code{vcov_model_based} are reporting-scale matrices; the corresponding
\code{vcov_optimizer_robust} and \code{vcov_optimizer_model_based} fields retain
the internal-coordinate versions.
Contemporaneous one- and two-step bootstraps transform every replicate before
calculating reported SEs and covariances.
}

\examples{
    mysummary <- simulate_data(n=1000, prep_asymptotic_se=FALSE)
    mymodel <- MCMmodel(mysummary, kernel="contemporaneous")
    my_result <- MCMfit(mymodel, mysummary, compute_se=FALSE,
                        optimizers="lbfgs", optim_iters=1, learning_rate=0.1)
    my_result_df <- as.data.frame(my_result)
    # Print loss
    print(my_result_df$loss)
}
