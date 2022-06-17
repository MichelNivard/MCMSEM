# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `torch-dev` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.  
This branch is for a potential move to torch for R backend

## Patch notes thus far (v0.3.1-dev-torch)
### Torch-specific (v0.3.1)
 - Moved to torch backend
 - Optimization is now performed using a combination of RPROP and LBFGS optimizers, instead of nlminb
 - Added MCMresultclass to hold MCMfit results
 - Can be converted to a dataframe using `as.data.frame(result)`, then it is in line again with that of the non-torch version
 - Note the torch implementation is slightly slower in the simplest use-case (2 variables + 1 confounder), but scales significantly better to more variables and confounders.
 - It is now possible to disable using skewness or kurtosis in larger models to improve perforamnce. Use either `use_skewness` and `use_kurtosis` arguments to `MCMfit()`.
### From v0.2.1-dev
 - Added mcmmodelclass, this class describes the layout of the MCMSEM model complete with parameter matrices, parameter/starting values, and bounds.
 - Addded MCMmodel wrapper function to enable easy creation of mcmmodelclass instances for users
 - Added MCMedit to make editing a model easier (e.g. adding or constraining parameters, changing bounds, etc.)
 - Added MCMfit
 - Setting the confounder to be negative is slightly harder than in the previous version, as it (for now) requires manual edits to the model:  
   `model <- MCMedit(model, "A", c(2,1), "-a1")`, note `c(2, 1)` are parameter coordinates in the A matrix and may depend on the number of confounders added to the model.
 - Output column names have changed to be identical to parameter names  
 - Added asymptotic calculation of standard errors for much faster runtimes

### Code updates
 - Update 17-06-2021 (labelled v0.3.1):
   - Removed several fixed TODO notes
     - Also removed notes about not being able to run MCMSEM for both positive and negative confounding, as this is not feasible in the current code format. It is easy enough now for users do this themselves through MCMedit.
   - Changed MCMedit argument names to more sensible ones
   - Moved `simulate_data` to MCMsimulate_data.R for more consistency in filenames.
   - Changed manual pages to new coding format
   - Changed test to new coding format
   - Added `use_skewness` and `use_kurtosis` arguments to `MCMfit()` to disable using one of them to estimate parameters in larger models.
   - Added `use_bounds` in MCMfit, if enabled, loss will be dramatically increased as parameter estimates get further away from bounds.
   - Added USAGE to readme
 - Update 15-06-2022-torch:
   - Fixed MCMfit so it now actually works
   - Merged changes to std.err from `dev`
   - For now, SE calculation is still done with R-matrices (similar to `dev`) as opposed to torch tensors, as I haven't found a way to make it work without significantly impacting performance.
   - Removed updates and todo that belong to `dev`
   - Added `silent` argument to `MCMfit()` to prevent printing loss at every step
   - Added MCMresult object to hold result as dataframe, loss, and history (all loss values, and model used)
   - Added `optim_iters` argument to `MCMfit()` to enable changing number of iterations of each optimizer

### Things still TODO:
1. Find a way to change .std.err and/or .jac.fn to use torch_tensors, if that's faster than R-matrices?

## Patch notes
- v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
- v0.1.0 - Initial commit

## Citation
If you use this package please include the following citation:
**#TODO: add reference**


## Installation

Currently this packge is not listed on CRAN and should therefore be installed from GitHub idrectly.
```
library(devtools)
install_github("https://github.com/matthijsz/MCMSEM")
```

## Usage

### Initial setup
After installing the package you can load it and start using its functionality. Note that the first time you load the package might take a while as torch installs several additional resources on first load.
``` 
library(MCMSEM)
```

Load your data into R, and prepare it for processing in MCMSEM, only functional variables should be passed to MCMSEM functions, so make a copy of your dataframe with variables like identifiers, timestamps and any character columns removed.
```
raw_data <- read.csv("myfile.csv", stringsAsFactors=FALSE)
my_data <- my_data[, -c("Identifier", "timestamp")] # This is just an example, column names will of course vary
```

### Creating a model
Next create an MCMmodel from your data
``` 
my_model <- MCMmodel(data)
```

By default, this will create model with one latent confounder, and constrained loadings of this confounder on your variables (i.e. 1 parameter for confounding). You can change this using the appropriate arguments.
``` 
my_extensive_model <- MCMmodel(data, n_confounding=2, constrained_a=FALSE)
```

By seting `constrained_a` to `FALSE` one parameter is fitted for each latent confounder for each variable (so n_confounding*n_variables confounders).

This model object contains the underlying matrices, parameter names and parameter values of the MCMSEM model, as well as some metadata. You can see the layout of the parameter matrices by running
```
print(my_model)
```

### Changing a model
You can edit several things about the model, such as freeing parameters, fixing parameters or starting values. Always use MCMedit for this purpose, **do not** edit the matrices by hand, as MCMedit parses your changes and ensures everything in the MCMmodel is changed accordingly. Below are some examples of what you may want to change.

In matrix A; separate a1 and a2, in the simplest case with two input variables this would be the same as setting `constrained_a` to `FALSE` in `MCMmodel`.
```
my_model <- MCMedit(my_model, "A", "a1", c("a1", "a2"))
```
In matrix Fm; freely estimate f1 variance, i.e. make it a parameter rather than a fixed value. Note the `c(1, 2)` in this function call refers to the coordinates of this parameter in the Fm matrix.
```
my_model <- MCMedit(my_model, "Fm", c(1, 2), "fm1")
```
Constrain b1_1 to zero
```
my_model <- MCMedit(my_model, "A", "b1_1", 0)
```

Change the starting value of b2 to 0.8
```
my_model <- MCMedit(my_model, "start", "b2", 0.8)
```
Set the starting value of all a parameters to 0.5
```
my_model <- MCMedit(my_model, "start", "a", 0.5)
```
Change the upper bound of a1 to 2
```
mcmmodel <- MCMedit(mcmmodel, "ubound", "a1", 2)
```
Set lower and upper bounds of a1 to c(-1, 1)
```
mcmmodel <- MCMedit(mcmmodel, "bound", "a1", c(-1, 1))
```
Set lower and upper bounds of all b parameters to -1, 1
```
mcmmodel <- MCMedit(mcmmodel, "bound", "b", c(-1, 1))
```
Set upper bound of all k parameters to 200
```
mcmmodel <- MCMedit(mcmmodel, "ubound", "k", 200)
```
Check the current bounds
```
mcmmodel$bounds
```
> **__NOTE__**: If you are using bounds, ensure `use_bounds` is set to TRUE in `MCMfit`
### Fitting a model
Now that you have applied the changes you want you can fit the model. You will have the pass your data to the `MCMfit` function again, we could store your data in the MCMmodel object, but that would result in unnecessary copies of potentially large datasets (not ideal).
```
my_result <- MCMfit(my_model, data)
```

The resulting `mcmresultclass` object contains your estimates and standard erros by default, as well as the final optimizer loss, loss history, and a copy of the model object that was used to generate the results.
Even though running `print(my_result)` produces output that looks like a dataframe, the object itself is not. If you want to continue working with your parameter and SE estimates, or write them to a file, you have to convert it to a dataframe first.
```
my_result_df <- as.data.frame(my_result)
```

By default, SEs are computed asymptotically, you can also compute the SEs through bootstrapping if you wish. Be warend, however, as this will take considerably longer.
``` 
my_result <- MCMfit(my_model, data, se_type='two-step', bootstrap_iter=100, bootstrap_chunks=1000)
```

If you are unsure the optimization is good enough, and wonder if a minimum is found, you can track the loss progress either by checking `my_result$history$loss` after the model completed, or by running
``` 
my_result <- MCMfit(my_model, data, silent=FALSE)
```

If needed you can change the number of iterations, and/or the optimizers learning rate. Note that two different optimizers are used sequentially. First RPROP, then LBFGS. Only the learning rate of RPROP can be changed.
``` 
my_result <- MCMfit(my_model, data, optim_iters=c(100, 12), learning_rate=0.03)
```

In large models using both skewness and kurtosis for parameter estimation may not be required, and disabling either of these can significantly improve performance. Note that this is still relatively untested though.
``` 
my_result_no_skew <- MCMfit(my_model, data, use_skewness=FALSE)
my_result_no_kurt <- MCMfit(my_model, data, use_kurtosis=FALSE)
```
