# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `dev` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.

Patch notes thus far (v0.2.0-dev):
 - Added mcmmodelclass, this class describes the layout of the MCMSEM model complete with parameter matrices, starting values, and bounds.
 - Addded MCMmodel wrapper function to enable easy creation of mcmmodelclass instances for users
 - Added mcmedit to make editing a model easier (e.g. adding or constraining parameters, changing bounds, etc.)
 - Update 16-05-2022: MCMfit now seems to work
   - TODO: Add progress bar (especially for bootstraps)
 - TODO: Test if results are identical to main branch
 - TODO: Enable setting confounder paths (a) to positive/negative

TODO: Move these semi-improvised notes to README and/or manual:
```
# Create model
mcmmodel <- MCMmodel(n_p = 2, n_f = 1, constrained_a=TRUE)  # n_p is n_phenotypes, will be replaced with data eventually
# This creates instance of reference class mcmmodelclass
# This instance contains:
#   - named matrices: just matrices with the labels -> b1_2 is effect of x1 on x2, b1_3 is effect of x1 on x3 etc.
#   - numeric matrices: starting values initially, but can be updated in optimization
#   - metadata: things like n_phenotypes, n_confounders etc, can hold as much as needed, as we'll probably need that stuff in other functions too
#   - param_values: current values of the parameters, so starting values initially but can be updated
#   - param_names: names of the parameters that can be changed
#   - param_coords: coordinates of param_values in their respective matrices (so the matrices can be easily updated during optimization)

# Can be edited via mcmedit, mcmedit in turn is designed such that it automatically updates the whole model
# Some examples:
# In matrix A; separate a1 and a2
mcmedit(mcmmodel, "A", "a1") <- c("a1", "a2")
mcmmodel # Check all matrices
# (In matrix Fm:) free f1 variance, i.e. make it a parameter rather than hard value
mcmedit(mcmmodel, "Fm", c(1, 2)) <- "fm1"
mcmmodel
# Constrain b1_1 to zero
mcmedit(mcmmodel, "A", "b1_1") <- 0
mcmmodel
mcmmodel$param_names  # The names of parameters that will be estimated (eventually)
# Now returns:
# [1] "a1"  "a2"  "b2"  "Fm1" "s1"  "s2"  "sk1" "sk2" "k1"  "k2"
# To change bounds...
# Check bounds (filled with defaults)
mcmmodel$bounds
# Change lower bound of a1 to -.5
mcmedit(mcmmodel, "lbound", "a1") <- -.5
mcmmodel$bounds
# Change the upper bound of 1 to 2
mcmedit(mcmmodel, "ubound", "a1") <- 2
mcmmodel$bounds
# Set lower and upper bounds of a1 to c(-1, 1)
mcmedit(mcmmodel, "bound", "a1") <- c(-1, 1)
mcmmodel$bounds
# Set lower and upper bounds of all b parameters to -1, 1
mcmedit(mcmmodel, "bound", "b") <- c(-1, 1)
mcmmodel$bounds
# Set upper bound of all k parameters to 200
mcmedit(mcmmodel, "ubound", "k") <- 200
mcmmodel$bounds
# a2 is added, Fm1 is added, and b1 is removed
# This doesn't work yet but the idea is that this is then passed to a fit function, e.g.
MCMfit(mcmmodel, estimate_SE=TRUE, n_iters=2000)
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

Use the `MCMSEM` function to fit MCMSEM models:
```
library(MCMSEM)
result <- MCMSEM(data)
```
Several options can be specified, such as the use of positive or negative confounding (or both):
```
library(MCMSEM)
result_positive <- MCMSEM(data, confounding='positive')
result_negative <- MCMSEM(data, confounding='negative')
results_combined <- MCMSEM(data, confounding='both')
```
**#TODO: add some text about what to do when**

By default Standard Errors for each estimate are computed through bootstrap, this can be disabled
```
result <- MCMSEM(data, compute_se=FALSE)
```

If computing SE takes to long, an alternative is to reduce the number of iterations (default is 200)
```
result <- MCMSEM(data, bootstrap_iter=100)
```

As a default a faster, but slightly less precise, two-step bootstrap is performed this behaviour can be changed if needed.
```
result <- MCMSEM(data, bootstrap_type='one-step')
```
**#TODO: maybe add some text about why this is faster and why this is (probably) just as accurate**

For testing purposes, this package also includes a function to simulate data which will work directly with the `MCMSEM` function:
```
data <- simulate_data(n=500000)
result <- MCMSEM(data)
```

### Patch notes
- v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
- v0.1.0 - Initial commit

