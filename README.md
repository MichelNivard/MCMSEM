# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## Citation
If you use this package please include the following citation:
**#TODO: add reference**


## Installation

Currently this packge is not listed on CRAN and should therefore be installed from GitHub directly.
```
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
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

