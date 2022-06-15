# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Development branch
Note this is the `torch-dev` branch, and **not** intended for end-users. If you would like to use MCMSEM yourself, please go to the main branch. If you would like to contribute to the code, feel free to check this branch out.  
This branch is for a potential move to torch for R backend

### Patch notes thus far (v0.2.0-dev-torch)
 - Moved to torch backend
 - Optimization is now performed using a combination of RPROP and LBFGS optimizers, instead of nlminb

### Code updates
 - Update 15-06-2022-torch:
   - Fixed MCMfit so it now actually works
   - Merged changes to std.err from `dev`
   - For now, SE calculation is still done with R-matrices (similar to `dev`) as opposed to torch tensors, as I haven't found a way to make it work without significantly impacting performance.
   - Removed updates and todo that belong to `dev`
   - Added `silent` argument to `MCMfit()` to prevent printing loss at every step

### Things still TODO:
1. Find a way to change .std.err and/or .jac.fn to use torch_tensors, if that's faster than R-matrices
