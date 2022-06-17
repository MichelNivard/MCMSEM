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

### Patch notes
- v0.1.1 
  - Added some TODO labels, added `'both'` option to `confounding`  argument which will run MCMSEM twice, once with negative, once with positive confounding and return both results. 
  - Added `bootstrap_chunk` argument. 
  - Added automatic standardization in `MCMSEM`
  - Added backups for internal `.m3m2v` and `.m4m2v` functions (should they ever be necessary).
- v0.1.0 - Initial commit