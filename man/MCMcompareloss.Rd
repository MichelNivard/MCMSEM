\name{MCMcompareloss}
\alias{MCMcompareloss}
\title{Compare MCM results}
\description{
  Compare loss of fitted MCM modelds given a specific dataset
}
\usage{
MCMcompareloss(results, test_data, weights=NULL, loss_type='auto',
               extensive_model_info=FALSE)
}
\arguments{
  \item{results}{list of MCM result objects. Names of list will be used as rownames of the output dataframe. If no names are provided they will be named model1, model2, ... etc}
  \item{test_data}{Either an MCM data summary object, or a matrix or dataframe with raw data, this test data / holdout data will be used to compute new loss}
  \item{weights}{If the models were fitted with weights, weights should be provided here as well.}
  \item{loss_type}{Type of loss to calculate, should be 'mse', 'smooth_l1' or c('mse','smooth_l1'). Alternatively 'auto' can be used in which case the loss type will be determined based on loss type used in the fitted models.}
  \item{extensive_model_info}{Returns more extensive details about the model, including other arguments such as device, low_memory, and also runtimes of the models.}
}
\value{
  A dataframe detailing the requested loss per model based on the provided data, as well as other statistics and details. Also includes loss difference column(s) detailing the difference of models 2-N with model1.
}

\examples{
    mydata <- simulate_data(n=3000, asdataframe=TRUE)

    # Split data in train and test sample
    split <- sample(c(1,2), nrow(mydata), prob=c(0.5, 0.5), replace=TRUE)
    mydata_train <- mydata[split == 1, ]
    mydata_test <- mydata[split == 2, ]

    # Generate and fit models on trainsample
    mydatasummary_train <- MCMdatasummary(
      mydata_train, prep_asymptotic_se=FALSE
    )
    mymodel <- MCMmodel(mydatasummary_train, kernel="contemporaneous")
    mymodel2 <- MCMmodel(mydatasummary_train, kernel="contemporaneous")
    my_result1 <- MCMfit(
      mymodel, mydatasummary_train, compute_se=FALSE,
      optimizers="lbfgs", optim_iters=1, learning_rate=0.1
    )
    my_result2 <- MCMfit(
      mymodel2, mydatasummary_train, compute_se=FALSE,
      optimizers="lbfgs", optim_iters=1, learning_rate=0.1
    )

    # compare loss in test sample
    MCMcompareloss(list(my_result1, my_result2), mydata_test)
}
