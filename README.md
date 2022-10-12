# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `dev-torch` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.  
As of version 0.4.0 it is possible to run MCMSEM on a GPU, see [MCMSEM on GPU](#mcmsem-on-gpu).

## Patch notes
### v0.16.0
 - Ported pre-calculation of `S.m` matrix in `MCMdatasummary()` to torch `jit_compile`, speeding up this part approximately 10x
   - ```
     nvars=30, nobs=    100   1000  10000
     S.m R-apply       1.4s  17.7s   155s
     S.m torch-jit     0.1s   1.1s  10.6s
     ```
   - Note this was one of the slowest part in `MCMdatasummary()` with many variables so I expect significant time savings
 - Moved all TorchScript code to list defined in `jit_funcs.R`
 - Added `jacobian_method` argument to `MCMfit` to change jacobian method for calculating standard errors
### v0.15.1
 - Moved definition of jit slownecker function to `.torch_fit` and `.std.err` to ensure compilation when the function is called, instead of when package is installed.
   - Thanks [@dfalbel](https://github.com/dfalbel) for your help.
### v0.15.0
 - Bugfix: Fixed improper calculation of K with diagonal S
 - Bugfix: set hardcoded values in K2 back to 3.0
### v0.14.1
 - Fixed improper assignment of `a` parameters, that would result in the left `A` matrix, instead of the right (see below)
   ```
    0  0  0  0  0  0  |  0  0  0  0  0  0
    0  0  0  0  0  0  |  0  0  0  0  0  0
    0  0  0  0  0  0  |  0  0  0  0  0  0
    a1 a2 a3 b...     |  a1 0  0  b...   
    a1 a2 a3 b...     |  a1 a2 0  b...   
    a1 a2 a3 b...     |  a1 a2 a3 b...   
    ```
### v0.14.0
 - Implemented `slowneckerproduct` when `low_memory=TRUE`, an iterative approach to `x %*% (y %x% y %x% y)`, which is significantly slower but does not require saving the entire `(y %x% y %x% y)` product, saving significant memory.
 - Added `monitor_grads` argument to `MCMfit()`. Setting `monitor_grads` to TRUE will make MCMSEM check for NaN gradients between each optimizer iteration, if any NaN gradients are found optimization is stopped early and a warning is printed. Then results are returned before the gradients are applied.
   - Currently `monitor_grads` does not work in LBFGS
   - Because this check will incur a performance penalty, `monitor_grads` is set to FALSE by default.
 - Moved cloud notes to `wiki/1.1 MCMSEM on the cloud.md`
 - Added wiki folder to `.Rbuildignore`
### v0.13.0
 - Implemented more efficient calculation of `K` from `S` when `S` is a diagonal matrix
   - Replaced `K %*% (S %x% S %x% S)` with `K * (diag(S) %x% diag(S) %x% (diag(S)))` when applicable 
### v0.12.1
 - Added `debug` argument to `MCMfit`, if enables prints detailed progress
### v0.12.0
 - Changed `MCMfit()` so it no longer uses `start_values` but uses `param_values`, instead. Now `start_values` are only stored for reference. This causes no change in behavior at initial `MCMfit()` call as through `MCMedit` `param_values` and `start_values` are always equal when `MCMfit()` is called for the first time.
 - Changed `MCMfit()` to also accept an `mcmresultclass` object to the `model` argument, so a user can train an already fitted model again (for more iterations, or lower learning_rate, etc.). 
   - ``` 
     data <- simulate_data()
     model <- MCMmodel(data)
     res <- MCMfit(model, data)
     res2 <- MCMfit(res, data)
     ```
### v0.11.0
 - Created `mcmstartvaluesclass` and changed all code for `MCMedit()`, `MCMmodel()` and `MCMfit()` accordingly
   - This is done so that `model$start_values["start", "a1_1"] <- 0` no longer works, as a measure to force users to change start values via `MCMedit("start", ...)`
 - Added option to edit all start values simultaneously to `MCMedit()`:
   - `newmodel <- MCMedit(model, "start", "all", new_values)`
   - Also added a brief description of this to wiki for `Editing an MCM model`
 - Added `outofbounds_penalty` argument to `MCMfit()` to increase penalty scaling when parameters are out bounds:
   - `loss *  2 ^ (distance_from_bounds * outofbounds_penalty)`
   - `MCMfit()` will produce an error if the value of `outofbounds_penalty` is set to a value < 0
   - If `outofbounds_penalty` is set to 0 `MCMfit()` will set `use_bounds` to `FALSE` instead, as that has the same effect but is more efficient
 - Added `outofbounds_penalty` to manual page for `MCMfit()`
 - Solved `unresolved reference n_obs` issue in `MCMcompareloss()
### v0.10.3
 - Fixed summary for models with 1 latent factor
 - Fixed summary for models without factor loadings and/or without causal paths 
### v0.10.2
 - Added `train_loss`, `train_chisq` and `train_bic` to `MCMcompareloss()` output, and renamed the newly calculated columns to `test_loss`, `test_chisq`, `test_bic`.
 - Added `train_n` to `MCMcompareloss(extensive_model_info=TRUE)`
   - Note that `train_n` can be different for result objects from prior versions. This is because previously the N was stored at `MCMmodel` not `MCMfit`. Now the N is stored at actual model fit `MCMfit()`.
 - Renamed `get_comoments()` to `.get_comoments()` to ensure it is a private function
 - Moved selection of loss function and optimizer function from the supported ones to their own functions in `local.R` named `.get_lossfunc(loss_type)` and `.get_optimizerfunc(optimizer)` respectively.
   - This is mainly done to allow for easier modification of the globally supported loss/optimizer lists
 - Completed transition from `for (i in 1:x)` to `for (i in seq_len(x))`
 - Allowed `use_kurtosis=FALSE` and `use_skewness=FALSE`
### v0.10.1
 - Renamed `data` argument to `test_data` in `MCMcompareloss()` to communicate the function is intended to be used with a separate holdout/test sample
 - Changed wording in manual page of `MCMcompareloss()` similarly
 - Changed example of `MCMcompareloss()` in the manual to reflect the intended use with test sample:
   - ```
     mydata <- simulate_data()
     # Split data in train and test sample
     split <- sample(c(1,2), nrow(data), prob=c(0.8, 0.2), replace=TRUE)
     mydata_train <- mydata[split == 1, ]
     mydata_test <- mydata[split == 2, ]
     
     # Generate and fit models on trainsample
     mymodel <- MCMmodel(mydata_train)
     mymodel2 <- MCMmodel(mydata_train)
     my_result1 <- MCMfit(mymodel, mydata_train)
     my_result2 <- MCMfit(mymodel2, mydata_train)
     
     # compare loss in test sample
     MCMcompareloss(list(my_result1, my_result2), mydata_test)
     ```
 - Added `.self$SE` to `mcmdatasummaryclass$copy()` 
### v0.10.0
 - Moved calculation of loss from `.torch_objective()` to separate `.calc_loss()` function as it is used in three different places
 - Added `weighted`, `loss_type`, and `optimizer` to `mcmresultclass$info`
 - Added `MCMcompareloss()`, a function to compare the loss (and other statistics) of different models, given a specific (holdout/test) dataset
   - > :warning: __WARNING__: Because of the addition of `weighted`, `loss_type`, and `optimizer` to `mcmresultclass$info` result objects created in earlier versions of MCMSEM (i.e. 0.9.1 or earlier) will not work in `MCMcompareloss()`
 - Added default starting values `(a=0.2, b=0, s=1, sk=0.8, k=4)` for new parameters created through `MCMedit()`, note it is still advised to set these yourself after creation of new paramters, but this way they are at least not all 0
 - Added manual pages for:
   - `MCMdatasumamry()`
   - `MCMsavesummary()`
   - `MCMcompareloss()`
   - `as.data.frame(mcmresultclass)`
   - `as.data.frame(mcmresultsummaryclass)`
   - `plot(mcmmodelclass)`
   - `plot(mcmresultclass)`
   - `print(mcmmodelclass)`
   - `summary(mcmresultclass)`
 - Updated manual pages for:
   - `MCMmodel()`
   - `MCMfit()`
 - Removed `Torch-specific` from versions in patch notes as the non-torch dev branch is deprecated and the torch-version is to be considered the only in-development branch
 - Merged `code updates` into `patch notes` for clarity
### v0.9.1
 - Fixed issue `object 'n' not found` in `MCMmodel()`
 - Fixed issue `object 'model' not found` in `MCMfit()`
### v0.9.0
 - Fixed an issue causing `summary(mcmresult)` to malfunction with causal paths between latents
 - Added `hdf5r` to NAMESPACE and Imports
 - Added `MCMdatasummary` to create an object with all required data for MCMSEM without storing full data. This could allow cohorts to more easily share this data without privacy concerns.
   - Note that this does make bootstrapping impossible
   - To ensure individual records cannot be retrieved from the summary the `S.m` matrix used for asymptotic SE calculation is processed into the required covariance matrix. Storing the `S.m` matrix at earlier stages would allow for faster data summary file generation, but could potentially allow for retrieval of individual records in some edge cases.
 - Moved calculation of `S.m` (in calculation of asymptotic SE) to `MCMdatasummary()` as it requires raw data.
 - Added option `prep_asymptotic_se` to `MCMdatasummary()` to allow for disabling pre-calculation of `S.m` matrix as it may result in unwanted compute time and storage space
 - Added `MCMsavesummary()` to save MCMdatasummaries to disk. File size estimates detailed below
 - Loading MCMdatasummaries is done as follows with the main summary function: `MCMdatasummary(path='mysummary.mcmdata')`
 - MCM data summaries stored to disk will have the suffix `.mcmdata` and are HDF5 files with the following structure:
   - m2 : second order co-moment matrix
   - m3 : third order co-moment matrix
   - m4 : fourth order co-moment matrix
   - meta/... : metadata, for now this contains `colnames`, `data_was_scaled`, `scale_data`, `n`, `ncol` and `colnames`, but this will automatically expand if we decide to add more meta_data to MCM data summary objects.
   - SE/computed : 1 if `prep_asymptotic_se` was set to true and `S.m` is pre-computed
   - sm : `S.m` matrix for asymptotic SE
   - SE/idx : list of different indices to select from the full `S.m` depending on `use_kurtosis` and `use_skewness`
 - Why HDF5 and not RData? Because HDF5...
   - Allows for more flexibility within R (i.e. `newobjectname <- MCMdatasummary(path=path)` instead of `load(path); newobjectname <- oldobjectname` where `oldobjectname` may not be known beforehand, or worse overlap with other objects)
   - Prevents storing class methods, so loading a file saved in an older version will still allow you to work with new methods (assuming that all required data is in the summary object)
   - Allows for easy and direct access from other software (though for now this is just future-proofing)
 - Moved warnings and errors related to the input dataframe from `MCMmodel` to `MCMdatasummary` as that is now used internally as well.
 - Changed internal structure of `MCMmodel`, `MCMfit` and underlying functions where necessary to also accept an `mcmdataclass` object as input instead of a dataframe
   - Simply put: if the input is a dataframe or matrix (i.e. raw data) `MCMmodel` and `MCMfit` will first generate an `MCMdatasummary` and continue working with that, to keep the code as simple as possible
 - Added `weights` argument to `MCMmodel`, `MCMdatasummary`, and `MCMfit` to allow for weighted MCMSEM (through weighted co-moment matrices)
   - Note it is recommended to make an `MCMdatasummary` first, then pass that object to `MCMmodel` and `MCMfit`, this will result in the co-moment-matrices being computed only ones, saving some time especially in larger datasets, and definitely when the data will be used in multiple different models.
 - Enabled setting `use_skewness` and `use_kurtosis` both to FALSE (with experimental warning).
 - Removed `Usage`, and `MCMSEM on GPU` text from `README.md`, instead referring to the wiki.
 -`.mcmdata` approximate file generation time and file size, with `prep_asymptotic_se` set to `TRUE`. Note these estimates are based on a single test run. Actual generation time will depend heavily on your system, and file size may vary depending on how well the summary can be compressed.  
Note that runtime is long, but that this is partly (or mostly, with many variables) due to pre-computation of the `S.m` matrix, so this upfront calculation will save time later on.

| N variables | generation time | file size |
|-------------|-----------------|-----------|
| 2           | 7.36 s          | 48.8 KB   | 
| 4           | 8.42 s          | 51.1 KB   |
| 6           | 10.76 s         | 70.1 KB   |
| 8           | 16.62 s         | 169 KB    |
| 10          | 27.10 s         | 538 KB    |
| 12          | 47.44 s         | 1.51 MB   |
| 14          | 82.44 s         | 4.05 MB   |
| 16          | 165.81 s        | 10.4 MB   |
| 18          | 278.62 s        | 22.9 MB   |
| 20          | 464.83 s        | 47.5 MB   |
| 22          | 730.05 s        | 91.8 MB   |
| 24          | 1222.58 s       | 170 MB    |
| 26          | 2108.23 s       | 295 MB    |
| 28          | 4072.52 s       | 505 MB    |

### v0.8.0
 - Fixed an issue that (depending on device) could cause fit to fail
 - Fixed `summ not found` issue in `summary(result)`
 - Changed optimizer iteration loops from `1:x` to `seq_len(x)` to allow for disabling one of the optimizers by setting its `optim_iters` to 0.
 - In `MCMfit` changed sqrt(S) to sign(S) * sqrt(abs(S)) to prevent NaN in cases of negative values (still in testing phase)
### v0.7.4
 - `loss` reported in result and summary objects is now the loss without bound scaling. You can still obtain the final loss with bound scaling via `res$history$loss[length(res$history$loss)]`
 - Added `...` to `plot(model, ..)` for additional arguments to be passed on to `qgraph`
 - Fixed a bug with `K` matrix when n_factors > 1 resulting from changes in v0.7.3
 - Fixed MCM summary so it now holds a deep copy of a result object, instead of a reference
 - Updated `wiki/4. MCMSEM results.md`
 - Added information on CPU threading to `1. Installing MCMSEM.md`
 - Updated manual pages based on recent changes
 - Added details on the layout of higher co-moment matrices to the `Advanced` section in `2.1 Creating an MCM model.md`.
### v0.7.3
 - Fixed an issue causing all bounds to be 0-100 (kurtosis-bounds) by default
 - Updated `wiki/3. Fitting an MCM model.md`
 - Values in `K` matrix for non-free kurtosis is now updated based on values in `S`
 - Allowed user to specify the first optimizer type through the `optimizer` argument, the following are included but all but RPROP are untested: `'rprop', 'sgd' ,'rmsprop', 'asgd', 'adam', 'adagrad', 'adadelta'`
 - Enabled models with 0 latent factors
 - Moved calculation of predicted M2, M3, M4 matrices to `.get_predicted_matrices` as the code was copy-pasted 3 times
### v0.7.2
 - Fixed summary(MCMresult) to also work when `use_skewness` or `use_kurtosis` is set to `FALSE`
 - Added `estimates` argument to `as.data.frame(MCMresultsummary)` to choose the estimates table to return (parameters, variances, skewness or kurtosis)
 - Added `wiki/2.0 Generate test data.md`
 - Updated `wiki/2.1 Creating an MCM model.md`
 - Updated `wiki/2.2 Editing an MCM model.md`
 - Added S3 method `print.mcmmodelclass` to allow for printing specific matrices, e.g. `print(mcmmodel, matrix="A")`
### v0.7.1
 - Fixed an issue that caused negative parameters to fail due to missing starting values
 - K matrix in `model$num_matrices$K` and `model$named_matrices$K` now properly display the value 3 as x,x,x,x-kurtosis for latent factors
 - Value of 3 for x,x,x,x-kurtosis of latent factors is no longer hardcoded in `MCMfit` but rather obtained from `model$num_matrices$K`
 - Added check for non-numeric columns and ID-column in `MCMmodel`
 - Changed stopping errors of multiple latent variables to warnings to allow for testing
 - Added observed and predicted comoment matrices to `result$history`
### v0.7.0
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
### v0.6.2
 - Added `wiki` folder with markdown documents to be used for the Wiki upon release.
 - Added argument `loss_type` to `MCMfit()` which allows users to change loss from MSE to smooth_L1. (regular L1 is also implemented but unlisted as it is untested)
 - Created `MCMSEMversion` variable in `local.R`, for easy access within the package (e.g. to store in result objects), and for easy access for users: `MCMSEM::MCMSEMversion`
 - Changed `License` field in `DESCRIPTION` to match CRAN standards
 - Replaced `for (i in 1:length(x))` with the safer `for (i in seq_along(x))` throughout
 - Replaced `for (i in 1:nrow(x))` with the safer `for (i in seq_len(nrow(x)))` throughout
 - Some basic cleanup without functional difference
### v0.6.1
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
### v0.6.0
 - Added `confounding_names` argument to `MCMmodel` to allow users to name latent factors
 - Original column names and latent factor names (user-provided or f1, f2, ..., fn) are now stored in `mcmmodel$meta_data`
 - Added `summary(mcmresult)` which returns a lavaan-partable-like dataframe using original column names and stored latent names
 - Added `plot(mcmresult)`which plots a qgraph visualization of the model using original column names and stored latent names
 - Added `plot(mcmmodel)` which allows visualizing the model before running `MCMfit()`, output is identical to `plot(mcmresult)` only with graph weights fixed to approximately 1.0
 - Changed check for pre-scaled data from checking `data == data_scaled` to checking if all columns mean and sd match expected values, significantly speeding up model generation with larger pre-scaled datasets
### v0.5.0
 - Changed call to R `cov()` in asymptotic SE calculation to custom torch solution (to be replaced with `torch_cov` once this is implemented), significantly improving performance.
 - Added `low_memory` option to `MCMfit`, which when enabled forces aggressive garbage collection during optimization. This can help run larger models on GPUs. Enabling `low_memory` does significantly impact performance, especially when using a CPU device, therefore it is off by default.
 - Added `info` field to `MCMresult` object which stores MCMSEM version, as well as arguments used to obtain the result.
### v0.4.3
 - Fixed `Error: invalid assignment for reference class field ‘param_values’`.
 - Added `int/CITATION` to allow users to run `citation("MCMSEM")`.
 - Renamed `MCMSEMmodelclass.R` to `MCMmodelclass` in line with other filenames.
 - Removed `.m2m2v`, `.m3m2v`, `.m4m2v`, `.torch_m2m2v`, `.torch_m3m2v` and `.torch_m4m2v` from `local.R` as they are no longer used (superseded by mask functions like `.torch_m2m2v_mask`)
 - Small memory tweaks in `MCMfit`
 - Moved to `torch nn_mse_loss` function significantly improving performance and reducing memory usage further
 - Added MCMedit option to change multiple parameters simultaneously using coordinates:
   - `mcmmodel <- MCMedit(mcmmodel, "A", list(6:10, 1), c("a11", "a12", "a13", "a14", "a15"))`
 - Added MCMedit option to remove all parameters of given type simultaneously:
   - `mcmmodel <- MCMedit(mcmmodel, "A", 'b', 0)`
### v0.4.2
 - Backend changes to standard error computation significantly improving performance.
 - Moved `model` object in `mcmresultclass` one level up so it can be accessed via `result$model`
 - Updated `model` object in `result` such that results are also stored in the matrices at appropriate locations
 - Changed `.jac.fn` to use torch tensors within the function significantly improving performance.
 - Changed `mcmresultclass` definition so `model` is now at top level of the object
 - Changed `mcmmodelclass` definition so an empty instance can be generated (for parsing/loading purposes)
 - The local `model$copy()` in `MCMfit` is now updated with resulting parameter values at the end of the function, so the returned model contains all parameters in the correct matrix.
 - Added [MCMSEM on GPU](#mcmsem-on-gpu) to README
### v0.4.1
 - Minor tweaks to `MCMfit` to slightly reduce (V)RAM usage
### v0.4.0
 - Added `device` argument to `MCMfit`, and reformated `MCMfit` to be device-agnostic. This enables using cuda (GPU)
 - Overhauled backend of `MCMfit` and `.objective` for improved performance:  
   - Optimizer-only runtime down from 30 minutes to 30 seconds (30 variables, 5 confounders, c(100, 25) iters, no SE, no kurtosis, no bounds)
 - Significant backend changes to `MCMfit` to significantly improve performance
   - (Nearly) everything ported to pure torch, no more indexing in `.objective`, no more `.m3m2v` loops etc.
   - Note this requires significantly more upfront work (creating base matrices, masks, etc), but results in ~60x performance boost in optimization
 - Enabled using with CUDA device
 - Fixed torch-implementation of quadratic scaling when parameters are out of bounds
 - Enabled using `use_bounds` with a CUDA device
 - Changed bootstrap code to be in line with current versions
### Code update 20-06-2021
 - Added `runtimes` element to MCMresult that gives some insight intor runtimes of several steps (mostly for development purpposes). 
 - Added `device` argument to `MCMfit` in preparation of CUDA-enabled version
### v0.3.1
 - Moved to torch backend
 - Optimization is now performed using a combination of RPROP and LBFGS optimizers, instead of nlminb
 - Added MCMresultclass to hold MCMfit results
 - Can be converted to a dataframe using `as.data.frame(result)`, then it is in line again with that of the non-torch version
 - Note the torch implementation is slightly slower in the simplest use-case (2 variables + 1 confounder), but scales significantly better to more variables and confounders.
 - It is now possible to disable using skewness or kurtosis in larger models to improve performance by reducing reliance on all co-moments. Set either `use_skewness` and `use_kurtosis` arguments to `MCMfit()` to `FALSE` to ignore specific co-moments.
 - Removed several fixed TODO notes
 - Also removed notes about not being able to run MCMSEM for both positive and negative confounding, as this is not feasible in the current code format. It is easy enough now for users do this themselves through MCMedit.
 - Changed MCMedit argument names to more sensible ones
 - Moved `simulate_data` to MCMsimulate_data.R for more consistency in filenames.
 - Changed manual pages to new coding format
 - Changed test to new coding format
 - Added `use_skewness` and `use_kurtosis` arguments to `MCMfit()` to disable using one of them to estimate parameters in larger models.
 - Added `use_bounds` in MCMfit, if enabled, loss will be dramatically increased as parameter estimates get further away from bounds.
 - Added USAGE to readme
### v0.3.0, Initial dev-torch
 - Fixed MCMfit so it now actually works
 - Merged changes to std.err from `dev`
 - For now, SE calculation is still done with R-matrices (similar to `dev`) as opposed to torch tensors, as I haven't found a way to make it work without significantly impacting performance.
 - Removed updates and todo that belong to `dev`
 - Added `silent` argument to `MCMfit()` to prevent printing loss at every step
 - Added MCMresult object to hold result as dataframe, loss, and history (all loss values, and model used)
 - Added `optim_iters` argument to `MCMfit()` to enable changing number of iterations of each optimizer
### v0.2.1-dev
 - Added mcmmodelclass, this class describes the layout of the MCMSEM model complete with parameter matrices, parameter/starting values, and bounds.
 - Addded MCMmodel wrapper function to enable easy creation of mcmmodelclass instances for users
 - Added MCMedit to make editing a model easier (e.g. adding or constraining parameters, changing bounds, etc.)
 - Added MCMfit
 - Setting the confounder to be negative is slightly harder than in the previous version, as it (for now) requires manual edits to the model:  
   `model <- MCMedit(model, "A", c(2,1), "-a1")`, note `c(2, 1)` are parameter coordinates in the A matrix and may depend on the number of confounders added to the model.
 - Output column names have changed to be identical to parameter names  
 - Added asymptotic calculation of standard errors for much faster runtimes
### v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
### v0.1.0 - Initial commit
 
### Things still TODO:
1. Figure out source of memory error with 30 variables and `low_memory=TRUE` 
2. Add gradient history to output
3. Move moment matrices in result to `res$observed` and `res$predicted`
4. Allow for use of multiple different optimizers: `optimizer=c("rprop", "rprop")`
5. Get `monitor_grads` to work in LBFGS
6. Create/update wiki/manual pages for:
   1. MCMdatasummary() - Include recommendation for generating datasummary object of data which will be used in different models
   2. MCMsavesummary()
   3. weighted analysis
   4. MCMcompareloss()
7. Torch-version of `S.m` generation
8. Add Hessian
9. Find a way to get the full jacobian using torch? (and have it be faster than the default jacobian with the current .jac.fn)
10. Expand checks in `MCMmodel` 

## Citation
If you use this package please include the following citation:  
Tamimy, Z., van Bergen, E., van der Zee, M. D., Dolan, C. V., & Nivard, M. G. (2022, June 30). Multi Co-Moment Structural Equation Models: Discovering Direction of Causality in the Presence of Confounding. [https://doi.org/10.31235/osf.io/ynam2](https://doi.org/10.31235/osf.io/ynam2)


## Installation

Currently this packge is not listed on CRAN and should therefore be installed from GitHub directly.
```
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
```

See the wiki `Instlaling MCMSEM` for more details.

## Usage

See the wiki, we recommend starting at `Generate test data`.

## MCMSEM on GPU

See the `Installing MCMSEM` wiki
