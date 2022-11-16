# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `dev-torch` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.  
As of version 0.4.0 it is possible to run MCMSEM on a GPU, see [MCMSEM on GPU](#mcmsem-on-gpu).

## Test observations
 - Anything above 25 LBFGS iterations tends to not matter at all
 - The minimum that is achieved depends far more on the LR and iterations of RPROP
   - Note however that the number of RPROP iterations quickly increases with larger models (1000-2000ish seemed to work well for 5-6 variables and 1 latent)
   - update: Scratch that, make it 10,000 for 7 variable models
 - Randomly generated a/b matrices used in `simulate_data` is unlikely to be reproduced by `MCMfit()`, my thoughts:
   - Random generation may be likely to generate completely unrealistic data
   - May generate data with multiple equally viable solutions
 - More variables leads to less accurate estimates (makes sense)
 - Variables with no relation to the others don't seem to have much of an impact (luckily)
 - Adding additional mutually exclusive latent factors (i.e. `a1_1` `a2_2`) to the 'true' model seems to make estimation significantly harder (especially for `a` parameters)
   - This may also result in unrealistic data though, I'm not entirely sure...
   - update: Seems likely, as once I increased the model to 7-8 variables these estimates became far more accurate. So my first tests with smaller datasets were just not identified.
 
## Patch notes
### v0.24.0
 - Fixed an issue causing the number of columns to be stored in results, instead of the number of samples (resulting in incorrect `BIC`)
 - Removed argument `silent` from `MCMfit()`, as its use was highly comparable to that of `debug`, loss printing during optimization is now included in `debug` and `silent` is removed
   - Additionally, enabled live loss-printing with CUDA devices (with the warning that it may significantly impact performance)
   - Updated `MCMfit()` manual page accordingly
 - Added wiki page `0.1 Rationale.md`
 - Updated wiki pages `Generating test data`, `Creating an MCM model`, `Editing an MCM model` and `MCMSEM results`, note these are currently still unfinished.
   - The main thing to do still is update the generated data and update the figures.
 - Added `debug` argument to `MCMsavesummary()`
   - Added `debug` entry to `MCMsavesummary()` manual page
 - Added wiki page 5 `Post hoc` with headers (unfinished), which will include:
   - MCM result plots
   - Loss history plots
   - Gradients/gradient plots
   - `MCMcompareloss()`
 - Moved patch-notes of version `0.21.1` and earlier to `PATCHNOTES.md` to clear readme file.
   - Note when this version is pushed to master all patch-notes will migrate there.
### v0.23.0
 - > :warning: __WARNING__: The default behavior of `simulate_data()` changed such that it now automatically returns an `mcmdatasummaryclass` object. To keep the raw simulated data in line with previous versions please use `simulate_data(asdtaframe=TRUE)`.
 - Added argument `asdataframe` to `simulate_data()`, and set the default to `FALSE`. this argument does the following:
   - `TRUE`: Returns the raw simulated data as data frame
   - `FALSE`: Returns an `mcmdatasummaryclass` object of the generated data
 - Added argument `...` to `simulate_data()`, used for passing other arguments to `MCMdatasummary()` when `asdataframe` is `FALSE`
 - Updated manual pages to be fully in line with v0.23.0 conventions for
   - `MCMcompareloss()`
   - `MCMdatasummary()`
   - `MCMedit()`
   - `MCMfit()`
   - `MCMmodel()`
   - `MCMsavesummary()`
   - `simulate_data()`
   - `as.data.frame(mcmresultclass)`
   - `as.data.frame(mcmresultsummaryclass)`
   - `plot(mcmmodelclass)`
   - `plot(mcmresultclass)`
   - `print(mcmmodelclass)`
   - `summary(mcmresultclass)`
### v0.22.0
 - Significantly changed `simulate_data` to allow for generation of N variables, in line with what the package is now capable of
   - Input `a1` changed to `a` which should now be an N_variables * N_latents matrix with `a` parameters
   - Inputs `b1` and `b2` changed to `b` which should now be an N_variables * N_variables matrix with `b` parameters and diagonal of 1
   - Added `skew` argument, boolean vector describing which variables should contain skewness
   - Added `kurt` argument, boolean vector describing which variables should contain kurtosis
   - `shape` and `df` still function the same, only are now used across all variables with skewness and kurtosis respectively
   - The default behavior of `simulate_data()` remains functionally the same
 - Added argument `parameterlegend` to `plot(result$gradients)` to allow for disabling legend, as this can sometimes obscure important details in the graph.
 - Changes in `MCMdatasummary`:
   - Changed possible ID column detection error to a warning
   - Changed `low_memory=0` behavior to use `covchunked` with 16 chunks as test results indicate this is faster with higher memory usage compared to base covariance
   - Changed the way data summary objects are saved from gzipped strings to H5 floats significantly speeding up the saving process (1.3s instead of 33s at 18 input variables), as well as saving space

 
### Things still TODO:
1. Test-code for (nearly) all configurations of MCMmodel/fit/etc
2. Create/update wiki pages for:
   1. MCMdatasummary() - Include recommendation for generating datasummary object of data which will be used in different models
   2. MCMsavesummary()
   3. weighted analysis
   4. MCMcompareloss()
   5. Gradient histories/`monitor_grads`
   6. Update manual `simulate_data()`
3. Add Hessian
4. Expand checks in `MCMmodel`
5. Get `monitor_grads` to work in LBFGS
6. Find a way to get the full jacobian and/or hessian using torch, and have it be faster than the default jacobian with the current .jac.fn
   - This is not feasible in the current implementation as for MCMSEM additional non-variable input arguments to `jacobian`/`hessian` are required (i.e. the fixed format matrices), other behavior-changing arguments like `low_memory` can be worked around by creating different functions for each type, but the matrices cannot be hardcoded. `torch.autograd.functional.jacobian`/`torch.autograd.functional.hessian` (not implemented in R torch but could theoretically be used via TorchScript) does not allow for non-variable (i.e. non-grad) arguments. The solution in Python is to use a class to hold fixed format objects, but custom classes are not available in TorchScript.
   - Note for future development, this will be possible if/when:
     1. R torch adds a more flexible version of `torch.autograd.functional.jacobian`/`torch.autograd.functional.hessian` which allows for additional non-variable input arguments, unlikely to happen any time soon though.
     2. TorchScript is expanded to allow for inclusion of custom classes. This seems even more unlikely than the previous point.
     3. MCMSEM migrates to a full Python backend (this is up to you)
     4. We create our own jacobian/hessian function in TorchScript from scratch

## Citation
If you use this package please include the following citation:  
Tamimy, Z., van Bergen, E., van der Zee, M. D., Dolan, C. V., & Nivard, M. G. (2022, June 30). Multi Co-Moment Structural Equation Models: Discovering Direction of Causality in the Presence of Confounding. [https://doi.org/10.31235/osf.io/ynam2](https://doi.org/10.31235/osf.io/ynam2)


## Installation

Currently, this package is not listed on CRAN and should therefore be installed from GitHub directly.
```
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
```

See the wiki `Instlaling MCMSEM` for more details.

## Usage

See the wiki, we recommend starting at `Generate test data`.

## MCMSEM on GPU

See the `Installing MCMSEM` wiki
