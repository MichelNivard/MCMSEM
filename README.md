# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## MCMSEM version 0.25.0
Welcome to the new and improved MCMSEM. If you want to use the MCMSEM version as it was used in [the original publication](https://doi.org/10.31235/osf.io/ynam2), please go to the [v0.1.1 release](https://github.com/zenabtamimy/MCMSEM/releases/tag/v0.1.1).

This version is considerably more powerful than our previous version. Some highlights:
 - Expanded to allow for any N variables instead of just two (use at your own risk)
 - Far more flexibility for custom model creation via MCMSEM model objects
 - More detailed fit statistics through custom MCMSEM result objects
 - Significantly improved performance, and enabled optimization on GPU
 - Asymptotic calculation of standard errors (bootstrapping no longer required)
 - Exportable data, making it easier for researchers to share moment matrices for MCMSEM without sharing raw data

If you are new to this version of MCMSEM we highly recommend reading our Wiki before starting, as the syntax for using MCMSEM has changed significantly since `v0.1.1`.

## Citation
If you use this package please include the following citation:  
Tamimy, Z., van Bergen, E., van der Zee, M. D., Dolan, C. V., & Nivard, M. G. (2022, June 30). Multi Co-Moment Structural Equation Models: Discovering Direction of Causality in the Presence of Confounding. [https://doi.org/10.31235/osf.io/ynam2](https://doi.org/10.31235/osf.io/ynam2)


## Installation

Currently, this package is not listed on CRAN and should therefore be installed from GitHub directly.
```
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
```

See the wiki `Installing MCMSEM` for more details.


### MCMSEM on GPU

See the `Installing MCMSEM` wiki

## Usage

See the wiki, we recommend starting at `Rationale`. In oversimplified terms, the general syntax is as follows:
```
data <- read.csv("mydata.csv")       # Read in your data
mcmdata <- MCMdatasummary(data, ...) # Convert to MCM data summary
base_mod <- MCMmodel(mcmdata, ...)   # Create an MCM model
mod <- MCMedit(mcmdata, ...)         # Edit an MCM model (if needed)
res <- MCMfit(mod, mcmdata, ...)     # Fit an MCM model
summary(res)                         # Generate output (if needed)
plot(res)                            # Generate visualizations (if needed)
```

## Contribute

If you would like to contribute to MCMSEM, please do so via the [dev-torch branch](https://github.com/zenabtamimy/MCMSEM/tree/dev-torch). 