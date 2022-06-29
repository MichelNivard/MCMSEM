# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## DEPRECATED development branch
Note this is the old `dev` branch of version 0.2 of MCMSEM, this is **not** intended for end-users.  
If you would like to use MCMSEM yourself, please go to the `main` branch.   
If you would like to contribute to the code, please go the `dev-torch` branch as development continues there.

### Patch notes thus far (v0.2.1-dev)
 - Added mcmmodelclass, this class describes the layout of the MCMSEM model complete with parameter matrices, parameter/starting values, and bounds.
 - Addded MCMmodel wrapper function to enable easy creation of mcmmodelclass instances for users
 - Added MCMedit to make editing a model easier (e.g. adding or constraining parameters, changing bounds, etc.)
 - Added MCMfit
 - Setting the confounder to be negative is slightly harder than in the previous version, as it (for now) requires manual edits to the model:  
   `model <- MCMedit(model, "A", c(2,1), "-a1")`, note `c(2, 1)` are parameter coordinates in the A matrix and may depend on the number of confounders added to the model.
 - Output column names have changed to be identical to parameter names  
 - Added asymptotic calculation of standard errors for much faster runtimes
 
### Code updates
 - Update 15-06-2022 (labelled v0.2.1):
   - Significantly improved performance of asymptotic SE calculation
   - Removed TODO "Add arguments for fitting either x->y or y->x path as opposed to both (which should remain the default)", this can easily be achieved through MCMedit, up to the user.
   - Changed MCMedit argument names to more sensible ones
   - Moved `simulate_data` to MCMsimulate_data.R for more consistency in filenames.
   - Changed manual pages to new coding format
   - Changed test to new coding format
 - Update 14-06-2022:
   - Added asymptotic computation of standard errors
   - Fixed bug with standardizing data 
 - Update 03-06-2022:
   - Fixed an issue causing error with dataframe input.
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
   - added `scale_data` option to `MCMmodel` to allow disabling automating data scaling
     - This option is stored in `model$meta_data$scale_data` and used again in MCMfit so the original dataframe is not changed
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
