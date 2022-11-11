\name{MCMfit}
\alias{MCMfit}
\title{Fit an MCM model}
\description{
  Fits an MCM sturctural equation model as defined in the MCMmodelclass object
}
\usage{
MCMfit(mcmmodel, data, weights=NULL, compute_se=TRUE, se_type='asymptotic', optimizers=c("rprop", "lbfgs"),
       optim_iters=c(50, 12), learning_rate=c(0.02, 1), loss_type='mse', bootstrap_iter=200, bootstrap_chunks=1000,
       silent=TRUE, use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE, device=NULL, device_se=NULL,
       low_memory=FALSE, outofbounds_penalty=1, monitor_grads=FALSE, debug=FALSE, jacobian_method='simple')
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
  \item{bootstrap_chunk}{Chunksize of two-step bootstrap. Lower numbers will lead to better performance but may lead to less accurate estimates.}
  \item{silent}{Run the analysis silently. If set to FALSE, will print the optimizer loss at every step}
  \item{use_bounds}{Use parameter bounds defined in the model, set this to FALSE in larger models to improve compute times}
  \item{use_skewness}{Use skewness and co-skewness to estimate parameters}
  \item{use_kurtosis}{Use kurtosis and co-kurtosis to estimate parameters}
  \item{device}{torch device object, can be used to run optimization on GPU; device=torch_device('cuda'), see "MCMSEM on GPU" in the README, defaults to using CPU.}
  \item{device_se}{torch device object used during calculation of SE (defaults to the same as device). Primarily used for changing device_se to a CPU device.}
  \item{low_memory}{Forces a more memory-friendly version of the optimization at the cost of runtime (note this cost gets VERY significant with large (>30 variables) models). Only set this option to TRUE when running into CUDA out of Memory errors.}
  \item{outofbounds_penalty}{loss penalty scaling when parameters are out of bounds: loss *  2 ^ (distance_from_bounds * outofbounds_penalty)}
  \item{monitor_grads}{Monitor the gradients during optimization and return gradients at each optimization. This will aditionally cause non-LBFGS optimization to stop early when NaN gradients are encountered, and returns the last non-NaN result instead.}
  \item{debug}{Print more process information}
  \item{jacobian_method}{Method used for the `jacobian()` function in asymptotic SE calculation}
}
\value{
  MCMresult object. This can be converted to a dataframe with as.data.frame(result)
}

\examples{
    mysummary <- simulate_data()
    mymodel <- MCMmodel(mysummary)
    my_result <- MCMfit(mymodel, mysummary)
    my_result_df <- as.data.frame(my_result)
    # Print loss
    print(my_result_df$loss)
}
