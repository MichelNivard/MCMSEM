\name{MCMfit}
\alias{MCMfit}
\title{Fit an MCM model}
\description{
  Fits an MCM sturctural equation model as defined in the MCMmodelclass object
}
\usage{
MCMfit(mcmmodel, data, compute_se=TRUE, se_type='asymptotic', optim_iters=c(50, 12), loss_type='mse', bootstrap_iter=200,bootstrap_chunks=1000,
                   learning_rate=0.02, silent=TRUE, use_bounds=TRUE, use_skewness=TRUE, use_kurtosis=TRUE, device=NULL, low_memory=FALSE)
}
\arguments{
  \item{model}{MCMmodelclass object}
  \item{data}{The matrix or dataframe that was used to generate the MCMmodelclass object}
  \item{compute_se}{Compute SE for the parameters}
  \item{se_type}{How the SE should be calculated. Should be one of 'asymptotic': ..., 'two-step':... or 'one-step':...}
  \item{optim_iters}{Number of iterations for RPROP and LBFGS optimizers. If convergence is not reached, start by increasing the number of iterations for RPROP.}
  \item{loss_type}{How to compute loss, should be one of 'mse' (mean squared errors), or 'smooth_l1' (smooth L1).}
  \item{bootstrap_iter}{Number of bootstrap iterations}
  \item{bootstrap_chunk}{Chunksize of two-step bootstrap. Lower numbers will lead to better performance but may lead to less accurate estimates.}
  \item{learning_rate}{learning rate of the RPROP optimizer. Increase if the optimizer has difficulty finding a minimum, decrease if loss quickly jumps to NA}
  \item{silent}{Run the analysis silently. If set to FALSE, will print the optimizer loss at every step}
  \item{use_bounds}{Use parameter bounds defined in the model, set this to FALSE in larger models to improve compute times}
  \item{use_skewness}{Use skewness and co-skewness to estimate parameters}
  \item{use_kurtosis}{Use kurtosis and co-kurtosis to estimate parameters}
  \item{device}{torch device object, can be used to run optimization on GPU; device=torch_device('cuda'), see "MCMSEM on GPU" in the README, defaults to using CPU.}
  \item{low_memory}{Forces aggressive garbage collection when enabled. Only set this option to TRUE when running into CUDA out of Memory errors.}
}
\value{
  MCMresult object. This can be converted to a dataframe with as.data.frame(result)
}

\examples{
    mydata <- simulate_data()
    mymodel <- MCMmodel(data)
    my_result <- MCMfit(mymodel, data)
    my_result_df <- as.data.frame(my_result)
    # Print loss
    print(my_result_df$loss)
}
