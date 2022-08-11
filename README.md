# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `dev-torch` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.  
As of version 0.4.0 it is possible to run MCMSEM on a GPU, see [MCMSEM on GPU](#mcmsem-on-gpu).

## Patch notes thus far (v0.7.3-dev-torch)
### Torch-specific (v0.7.3)
 - Fixed an issue causing all bounds to be 0-100 (kurtosis-bounds) by default
 - Updated `wiki/3. Fitting an MCM model.md`
 - Values in `K` matrix for non-free kurtosis is now updated based on values in `S`
 - Allowed user to specify the first optimizer type through the `optimizer` argument, the following are included but all but RPROP are untested: `'rprop', 'sgd' ,'rmsprop', 'asgd', 'adam', 'adagrad', 'adadelta'`
 - Enabled models with 0 latent factors
### Torch-specific (v0.7.2)
 - Fixed summary(MCMresult) to also work when `use_skewness` or `use_kurtosis` is set to `FALSE`
 - Added `estimates` argument to `as.data.frame(MCMresultsummary)` to choose the estimates table to return (parameters, variances, skewness or kurtosis)
 - Added `wiki/2.0 Generate test data.md`
 - Updated `wiki/2.1 Creating an MCM model.md`
 - Updated `wiki/2.2 Editing an MCM model.md`
 - Added S3 method `print.mcmmodelclass` to allow for printing specific matrices, e.g. `print(mcmmodel, matrix="A")`
### Torch-specific (v0.7.1)
 - Fixed an issue that caused negative parameters to fail due to missing starting values
 - K matrix in `model$num_matrices$K` and `model$named_matrices$K` now properly display the value 3 as x,x,x,x-kurtosis for latent factors
 - Value of 3 for x,x,x,x-kurtosis of latent factors is no longer hardcoded in `MCMfit` but rather obtained from `model$num_matrices$K`
 - Added check for non-numeric columns and ID-column in `MCMmodel`
 - Changed stopping errors of multiple latent variables to warnings to allow for testing
 - Added observed and predicted comoment matrices to `result$history`
### Torch-specific (v0.7.0)
> :warning: __WARNING__: The argument names `n_confounding` and `confounding_names` in `MCMmodel()` have been changed to `n_latent` and `latent_names` in this version. Please update your code accordingly.
 - Changed references to latent factors from `confounding` to `latent` script-wide, for consistent labelling. 
 - Added 8 new arguments to `MCMmodel()` to change broad model specifications directly. Note for now these do not alter the initial model generation, but use `MCMedit()` after the default model is generated, slightly slowing down `MCMmodel`, so if we get to the point where giant models (i.e. >50 variables) are possible, this should be changed.
   - causal_observed: adds free parameters for causal paths between observed variables (default TRUE)
   - var_observed: adds free parameters for variance of observed variables (default TRUE)
   - skew_observed: adds free parameters for skewness of observed variables (default TRUE)
   - kurt_observed: adds free parameters for kurtosis of observed variables (default TRUE)
   - causal_latent: adds free parameters for causal paths between latent factors (default FALSE)
   - var_latent: adds free parameters for variance of latent factors (default FALSE)
   - skew_latent: adds free parameters for skewness of latent factors (default FALSE)
   - kurt_latent: adds free parameters for kurtosis of latent factors (default FALSE)
 - Changed `MCMedit()` to accept multiple parameter names, e.g.: `MCMedit(model, "A", c("b1_2", "b2_1"), 0)`. Note that all these parameters should still originate from the same matrix
### Torch-specific (v0.6.2)
 - Added `wiki` folder with markdown documents to be used for the Wiki upon release.
 - Added argument `loss_type` to `MCMfit()` which allows users to change loss from MSE to smooth_L1. (regular L1 is also implemented but unlisted as it is untested)
 - Created `MCMSEMversion` variable in `local.R`, for easy access within the package (e.g. to store in result objects), and for easy access for users: `MCMSEM::MCMSEMversion`
 - Changed `License` field in `DESCRIPTION` to match CRAN standards
 - Replaced `for (i in 1:length(x))` with the safer `for (i in seq_along(x))` throughout
 - Replaced `for (i in 1:nrow(x))` with the safer `for (i in seq_len(nrow(x)))` throughout
 - Some basic cleanup without functional difference
### Torch-specific (v0.6.1)
 - Changed `cat` calls in non-silent operation of `MCMfit` to use CR (`\r`) instead of newline (`\n`) to prevent flooding the console 
 - Split `local.R` into `local_gen_matrices.R`, `local_stderr.R`, and `local_torch_matrices.R`, also moved fit and objective functions from `MCMfit.R` to `local_fit.R`. Each script now contains local functions that are (mainly) used for their respective function names
 - Removed `knot` column from `summary(mcmresult)` as it is currently unused
 - Removed `par` column from `summary(mcmresult)` as it did not provide additional information
 - Changed `edge` column for skewness and kurtosis parameters to `~~~` and `~~~~` respectively
 - Renamed `MCMresultclass.R` to  `MCMresultclasses.R` as it now also holds `mcmresultsummaryclass` (putting it in a separate script didn't work)
 - Added `mcmresultsummaryclass`, a custom class now returned from `summary(mcmresult)` in order to store additional information. The `mcmresultsummaryclass` contains the following:
   - `df`: dataframe with parameter table (this used to be the only thing returned from `summary(mcmresult)`)
   - `result`: A reference to the result object that the summary is generated from
 - `as.data.frame(summary(mcmresult))` simply extracts the parameter table
 - `plot(summary(mcmresult))` is the same as `plot(mcmresult)`, but it is implemented as I suspect people will expect it to work
 - Added `n_obs` to `mcmmodel$meta_data`
 - Added several parameters to MCM result summary, most notably `chisq`, `aic` and `bic`
 - Changed `n_par` calculation in `summary(result)` to subtract parameters in skewness/kurtosis matrix if either of those was disabled during `MCMfit`
 - switched to Huber loss for stability (`nn_smooth_l1_loss()`), TODO: Make an optional switch between MSE and Huber loss
### Torch-specific (v0.6.0)
 - Added `confounding_names` argument to `MCMmodel` to allow users to name latent factors
 - Original column names and latent factor names (user-provided or f1, f2, ..., fn) are now stored in `mcmmodel$meta_data`
 - Added `summary(mcmresult)` which returns a lavaan-partable-like dataframe using original column names and stored latent names
 - Added `plot(mcmresult)`which plots a qgraph visualization of the model using original column names and stored latent names
 - Added `plot(mcmmodel)` which allows visualizing the model before running `MCMfit()`, output is identical to `plot(mcmresult)` only with graph weights fixed to approximately 1.0
 - Changed check for pre-scaled data from checking `data == data_scaled` to checking if all columns mean and sd match expected values, significantly speeding up model generation with larger pre-scaled datasets
### Torch-specific (v0.5.0)
 - Changed call to R `cov()` in asymptotic SE calculation to custom torch solution (to be replaced with `torch_cov` once this is implemented), significantly improving performance.
 - Added `low_memory` option to `MCMfit`, which when enabled forces aggressive garbage collection during optimization. This can help run larger models on GPUs. Enabling `low_memory` does significantly impact performance, especially when using a CPU device, therefore it is off by default.
 - Added `info` field to `MCMresult` object which stores MCMSEM version, as well as arguments used to obtain the result.
### Torch-specific (v0.4.3)
 - Small memory tweaks in `MCMfit`
 - Moved to `torch nn_mse_loss` function significantly improving performance and reducing memory usage further
 - Added MCMedit option to change multiple parameters simultaneously using coordinates:
   - `mcmmodel <- MCMedit(mcmmodel, "A", list(6:10, 1), c("a11", "a12", "a13", "a14", "a15"))`
 - Added MCMedit option to remove all parameters of given type simultaneously:
   - `mcmmodel <- MCMedit(mcmmodel, "A", 'b', 0)`
### Torch-specific (v0.4.2)
 - Backend changes to standard error computation significantly improving performance.
 - Moved `model` object in `mcmresultclass` one level up so it can be accessed via `result$model`
 - Updated `model` object in `result` such that results are also stored in the matrices at appropriate locations
### Torch-specific (v0.4.1)
 - Minor tweaks to `MCMfit` to slightly reduce (V)RAM usage
### Torch-specific (v0.4.0)
 - Added `device` argument to `MCMfit`, and reformated `MCMfit` to be device-agnostic. This enables using cuda (GPU)
 - Overhauled backend of `MCMfit` and `.objective` for improved performance:  
   - Optimizer-only runtime down from 30 minutes to 30 seconds (30 variables, 5 confounders, c(100, 25) iters, no SE, no kurtosis, no bounds)
### Torch-specific (v0.3.1)
 - Moved to torch backend
 - Optimization is now performed using a combination of RPROP and LBFGS optimizers, instead of nlminb
 - Added MCMresultclass to hold MCMfit results
 - Can be converted to a dataframe using `as.data.frame(result)`, then it is in line again with that of the non-torch version
 - Note the torch implementation is slightly slower in the simplest use-case (2 variables + 1 confounder), but scales significantly better to more variables and confounders.
 - It is now possible to disable using skewness or kurtosis in larger models to improve performance by reducing reliance on all co-moments. Set either `use_skewness` and `use_kurtosis` arguments to `MCMfit()` to `FALSE` to ignore specific co-moments.
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
 - Update 01-07-2021:
   - Fixed `Error: invalid assignment for reference class field ‘param_values’`.
   - Added `int/CITATION` to allow users to run `citation("MCMSEM")`.
   - Renamed `MCMSEMmodelclass.R` to `MCMmodelclass` in line with other filenames.
   - Removed `.m2m2v`, `.m3m2v`, `.m4m2v`, `.torch_m2m2v`, `.torch_m3m2v` and `.torch_m4m2v` from `local.R` as they are no longer used (superseded by mask functions like `.torch_m2m2v_mask`)
 - Update 30-06-2021 (labelled v0.4.2):
   - Changed `.jac.fn` to use torch tensors within the function significantly improving performance.
   - Changed `mcmresultclass` definition so `model` is now at top level of the object
   - Changed `mcmmodelclass` definition so an empty instance can be generated (for parsing/loading purposes)
   - The local `model$copy()` in `MCMfit` is now updated with resulting parameter values at the end of the function, so the returned model contains all parameters in the correct matrix.
   - Added [MCMSEM on GPU](#mcmsem-on-gpu) to README
 - Update 29-06-2021 (labelled v0.4.1):
   - Minor RAM tweaks in `MCMfit`, storing fewer R objects saves some VRAM (particularly for CUDA devices).
 - Update 24-06-2021 (labelled v0.4.0):
   - Significant backend changes to `MCMfit` to significantly improve performance
     - (Nearly) everything ported to pure torch, no more indexing in `.objective`, no more `.m3m2v` loops etc.
     - Note this requires significantly more upfront work (creating base matrices, masks, etc), but results in ~60x performance boost in optimization
   - Enabled using with CUDA device
   - Fixed torch-implementation of quadratic scaling when parameters are out of bounds
   - Enabled using `use_bounds` with a CUDA device
   - Changed bootstrap code to be in line with current versions
 - Update 20-06-2021:
   - Added `runtimes` element to MCMresult that gives some insight intor runtimes of several steps (mostly for development purpposes). 
   - Added `device` argument to `MCMfit` in preparation of CUDA-enabled version
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
1. Find a way to get the full jacobian using torch? (and have it be faster than the default jacobian with the current .jac.fn)
2. Expand checks in `MCMmodel` 

## Patch notes
- v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
- v0.1.0 - Initial commit

## Citation
If you use this package please include the following citation:  
Tamimy, Z., van Bergen, E., van der Zee, M. D., Dolan, C. V., & Nivard, M. G. (2022, June 30). Multi Co-Moment Structural Equation Models: Discovering Direction of Causality in the Presence of Confounding. [https://doi.org/10.31235/osf.io/ynam2](https://doi.org/10.31235/osf.io/ynam2)


## Installation

Currently this packge is not listed on CRAN and should therefore be installed from GitHub idrectly.
```
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
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


## MCMSEM on GPU
> **__Preface__**: Using a GPU is only advantageous under certain circumstances. Generally, if your intended model contains 5 or fewer input variables, or if you intend to run a model without kurtosis it will likely not be worth it, though this will all depend on your model, data, CPU, GPU etc.

### Requirements
1. An NVIDIA CUDA-enabled GPU, see [NVIDIA's website](https://developer.nvidia.com/cuda-gpus).
   - On top of this you will likely need a GPU with a significant VRAM capacity. For reference, our tests are performed on a GTX1080Ti (11GB VRAM) which is limited to running skew+kurtosis models with 18-20 input variables and 5 confounders.
2. Installed NVIDIA CUDA toolkit version 11.3 (note it has to be **exactly** version 11.3)
    - Make sure that environmental variables CUDA_HOME and CUDA_PATH are set properly
3. Installed cuDNN 8.4 for CUDA toolkit 11.3 see [NVIDIA's website](https://docs.nvidia.com/deeplearning/cudnn/archives/cudnn-840/install-guide/index.html) (again, it has to be **exactly** version 8.4).
4. Installed torch after all the above conditions are met.
   - If you already have torch for R installed, remove it: `remove.packages("torch")`, restart R, and install it again.

### Usage
Before running MCMSEM, verify that CUDA is available to your torch installation. Note that the first time you run `library(torch)` might take a few minutes as torch will then download and install required gpu-torch libraries.
``` 
library(torch)
cuda_is_available()  # Should return TRUE
```

If CUDA is available, setup a CUDA device and pass this to MCMfit.
``` 
cuda_device <- torch_device("cuda")
res <- MCMfit(my_model, data, device=cuda_device)
```

Running `MCMfit` on a GPU can result in very significant runtime improvements, for instance from our testing a model with 15 input variables and 5 confounders, using both skewness and kurtosis took an hour on CPU compared to 20 minutes on GPU.  
This, however, is contingent upon enough available VRAM. If you get a `CUDA out of memory` error your VRAM is insufficient given your model, try using `MCMfit(..., low_memory=TRUE)`.  
We do our best to ensure the code uses as little VRAM as possible, but due to the nature of GPU computing VRAM requirements will likely remain much stricter than RAM requirements.
