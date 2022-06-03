# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `dev` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.

### Patch notes thus far (v0.2.0-dev)
 - Added mcmmodelclass, this class describes the layout of the MCMSEM model complete with parameter matrices, parameter/starting values, and bounds.
 - Addded MCMmodel wrapper function to enable easy creation of mcmmodelclass instances for users
 - Added MCMedit to make editing a model easier (e.g. adding or constraining parameters, changing bounds, etc.)
 - Added MCMfit
 - Setting the confounder to be negative is slightly harder than in the previous version, as it (for now) requires manual edits to the model:  
   `model <- MCMedit(model, "A", c(2,1), "-a1")`, note `c(2, 1)` are parameter coordinates in the A matrix and may depend on the number of confounders added to the model.
 - Output column names have changed to be identical to parameter names  

### Code updates
 - Update 03-06-2022:
   - Fixed an issue causing starting values to not be properly assigned to a-parameters with >1 confounding
   - Updated bootstrap to work with >2 input data columns
   - Updated bootstrap to work with >1 confounding
   - Fixed todo related to >1 confounding in `MCMmodel.R`:
     - `if n_confounding > 1: make sure f2->x1 is set to 0, i.e. f1 -> x1,x2,x3,x4, f2 -> x2,x3,x4, etc.`
   - All of MCMfit should now work with an arbitrary number of input data columns and confounders
   - Changed error in MCMfit when provided with 3 columns of data to a warning that this is still experimental
   - Changed column names of MCMfit output to match parameter names
   - Added start_values as an element of MCMmodel class, start values can now be changed via `MCMedit`:
     - `MCMedit(model, "start", "a1", 0.3)`
 - Update 31-05-2022:
   - MCMfit works
   - Added basic progress bar to MCMfit bootstrap
   - Changed `mcmedit` to `MCMedit` to be in line with capitalization of other functions
   - Changed usage of `MCMedit` from `MCMedit(model, matrix, cell) <- value` to `model <- MCMedit(model, matrix, cell, value)`
   - `model <- MCMedit(model, matrix, cell, value)` now also creates a copy of model instead of modifying inplace.
   - Removed old `MCMSEM` and `.fn` functions from code
   - Removed old Usage text from README
   - Added options for negative versions of parameters (mainly for `-a1` as negative confounder).
     - This is done through MCMedit for now: `model <- MCMedit(model, "A", c(2,1), "-a1")`, model parser will recognize paramters labelled with "-".
     - Testing this resulted in a near-0 estimate of a... Should verify if this is the same as master branch
     - In a future version we may want to write a wrapper function to do `model <- MCMedit(model, "A", c(2,1), "-a1")` in a more user-friendly way?

### Things still TODO:
1. Check if results of positive/negative confounder are identical to master branch
2. Rename columns of MCMfit output to better reflect column names of input data, e.g. `b_age_pheno` instead of `b1_2`
   - Easiest is probably to store these names somewhere and simply replace them.
3. Create detailed manual page for MCMedit
4. Change MCMedit argument names to sensible ones
5. Make bootstrap MCMfit run in parallel
   - Note: current progress bar is not suited for parallel bootstrap
6. Update tests to reflect new coding style.
7. Move these semi-improvised notes to README and/or manual:
```
# Create model
mcmmodel <- MCMmodel(data, n_confounding = 1, constrained_a=TRUE)  # n_p is n_phenotypes, will be replaced with data eventually
# This creates instance of reference class mcmmodelclass
# This instance contains:
#   - named matrices: just matrices with the labels -> b1_2 is effect of x1 on x2, b1_3 is effect of x1 on x3 etc.
#   - numeric matrices: starting values initially, but are updated during optimization
#   - metadata: things like n_phenotypes, n_confounders etc, can hold as much as needed, as we'll probably need that stuff in other functions too
#   - param_values: current values of the parameters, so starting values initially but can be updated
#   - param_names: names of the parameters that can be changed
#   - param_coords: coordinates of param_values in their respective matrices (so the matrices can be easily updated during optimization), as well as multipliers to account for user-negative parameters

# Can be edited via MCMedit, MCMedit in turn is designed such that it automatically updates the whole model
# Some examples:
# In matrix A; separate a1 and a2
mcmmodel <- MCMedit(mcmmodel, "A", "a1", c("a1", "a2"))
mcmmodel # Check all matrices
# (In matrix Fm:) free f1 variance, i.e. make it a parameter rather than hard value
mcmmodel <- MCMedit(mcmmodel, "Fm", c(1, 2), "fm1")
mcmmodel
# Constrain b1_1 to zero
mcmmodel <- MCMedit(mcmmodel, "A", "b1_1", 0)
mcmmodel
mcmmodel$param_names  # The names of parameters that will be estimated (eventually)
# Now returns:
# [1] "a1"  "a2"  "b2"  "Fm1" "s1"  "s2"  "sk1" "sk2" "k1"  "k2"
# To change bounds...
# Check bounds (filled with defaults)
mcmmodel$bounds
# Change lower bound of a1 to -.5
mcmmodel <- MCMedit(mcmmodel, "lbound", "a1", -.5)
mcmmodel$bounds
# Change the upper bound of 1 to 2
mcmmodel <- MCMedit(mcmmodel, "ubound", "a1", 2)
mcmmodel$bounds
# Set lower and upper bounds of a1 to c(-1, 1)
mcmmodel <- MCMedit(mcmmodel, "bound", "a1", c(-1, 1))
mcmmodel$bounds
# Set lower and upper bounds of all b parameters to -1, 1
mcmmodel <- MCMedit(mcmmodel, "bound", "b", c(-1, 1))
mcmmodel$bounds
# Set upper bound of all k parameters to 200
mcmmodel <- MCMedit(mcmmodel, "ubound", "k", 200)
mcmmodel$bounds
# Set the starting value of b2 to 1.5
mcmmodel <- MCMedit(mcmmodel, "start", "b2", 1.5)
# Set the starting value of all a parameters to 1.0
mcmmodel <- MCMedit(mcmmodel, "start", "a", 1.0)

# Fit the model
MCMfit(mcmmodel, estimate_SE=TRUE, bootstrap_iter=2000)
```

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

TODO

### Patch notes
- v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
- v0.1.0 - Initial commit

