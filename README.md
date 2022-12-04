# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## MCMSEM version 0.26.1
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

Below you will find a short rationale with usage examples of MCMSEM. For more detailed descriptions please visit the wiki.

### A rationale for MCMSEM: this is BIG, you should care!

Let's get going, in this very short pre-tutorial I'll convince you why you should read the entire tutorial, the paper(s) and consider MCMSEM for your projects. 
This is an advertorial, not a full review of the method with all its good and bad, that's left for the paper and the rest of this wiki.

```{r}
library(devtools)
install_github("https://github.com/zenabtamimy/MCMSEM")
library(MCMSEM)
library(lavaan)
```

The basic premise is that by modeling higher order co-moments, not just covariance MCMSEM can do incredible things. 
It can (for example) estimate directional causal effects in the presence of other causal effects in the opposite direction based on continuous variables collected in a cross-sectional and observational setting.
Given certain assumptions about the confusers hold even causal effects in the presence of confounding. So to proof that to you, let me simulate data from a dense network with bidirectional causal relations between variables. 
In the simulation we only use direct paths between variables to induce correlations, no latent variables are present.

```
b <- matrix(c(   1,    .3,   0, .25  ,.25, 
                 .05,   1,   0,  .3,  .25,
                 .35,  .4,   1,  .3,   .2,
                 .15, .05,  .1,   1,  .48,
                 0.2,0.25, 0.25,  0,    1), 5,5,byrow=T)
# Latent variables don't load on the indicators:
a <- matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 0, 0), ncol=2)

# use the MCMSEM internal simuation tool:
simmdata<- simulate_data(n=25000,a=a,b=b,shape=c(7, 0, 3, 4, 5), df=c(0, 8, 0, 10, 12),asdataframe = T)

#simulate holdout data from the exact same process!
simdata.holdout <-  simulate_data(n=10000,a=a,b=b,shape=c(7, 0, 3, 4, 5), df=c(0, 8, 0, 10, 12),asdataframe = T)
cor(simmdata)
```

```
              [,1]      [,2]      [,3]      [,4]      [,5]
    [1,] 1.0000000 0.4684255 0.5023147 0.4906726 0.4480491
    [2,] 0.4684255 1.0000000 0.4912504 0.4994899 0.4563662
    [3,] 0.5023147 0.4912504 1.0000000 0.4886721 0.4533532
    [4,] 0.4906726 0.4994899 0.4886721 1.0000000 0.4431923
    [5,] 0.4480491 0.4563662 0.4533532 0.4431923 1.0000000
```

So we generated network data, which gives rise to 5 correlated variables, we simulated 25000 observations of 5 variables, we made sure these variables have some skewness and kurtosis. 
Then we computed the correlations between the variables, these correlations seem sort of consistent with the influence of a single latent variable (all 5 variables are correlated about equally). 
This is a known problem right? the data don't really identify a specific model, and might fit the wrong model rather well, actually. If we fit a latent variable model to these data in `lavaan` what happens?
```
### Lavaan naiveness:
# specify a single factor model:
model <- " F1=~ V1 + V2+ V3+ V4+ V5" 

#Fit a single factor modle to the data:
single.factor.model <- sem(model,data = simmdata)

### Based on respectable fit indices, the single factor model has pretty good fit....
fitmeasures(single.factor.model,fit.measures = c("cfi","rmsea"))
```

```
      cfi rmsea 
    0.998 0.022 
```

So that's a pretty solid fit to 25000 data points which aren't normally distributed! I think may people would happily accept the (wrong) model fits the data well in this case.
Let's look at the estimated parameters, I took liberty of omitting a part of the result (here and in other examples below) to improve readability of the page:

```
summary(single.factor.model,standardize=T)
```

```
    lavaan 0.6-12 ended normally after 23 iterations

      Estimator                                         ML


    Latent Variables:
                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
      F1 =~                                                                 
        V1                1.000                               0.801    0.695
        V2                0.998    0.011   92.422    0.000    0.799    0.697
        V3                1.098    0.012   93.495    0.000    0.879    0.707
        V4                1.391    0.015   92.872    0.000    1.114    0.701
        V5                1.183    0.014   86.523    0.000    0.947    0.643

    Variances:
                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
       .V1                0.687    0.008   88.249    0.000    0.687    0.517
       .V2                0.676    0.008   87.957    0.000    0.676    0.514
       .V3                0.771    0.009   86.480    0.000    0.771    0.499
       .V4                1.281    0.015   87.356    0.000    1.281    0.508
       .V5                1.273    0.014   94.116    0.000    1.273    0.587
        F1                0.641    0.011   57.069    0.000    1.000    1.000
```

### Let's try MCMSEM!
Okay we can fit the exact same model in MCMSEM, this will take longer, but that's because we are using MCMSEM for a simple model while it's meant for way more complex models...

```
# Prepare the data:
simmdatasumm <- MCMdatasummary(simmdata)

# specify the MCMSEM single factor model:
mod.single.fac <- MCMmodel(simmdatasumm, n_latent=1,
                          causal_observed = F ,constrained_a = FALSE)
res.single.fac  <- MCMfit(mod.single.fac , simmdatasumm,
                          optimizers=c("rprop", "lbfgs"), optim_iters=c(2000, 25), 
                          learning_rate=c(0.15, 1.0), monitor_grads = TRUE)
summary(res.single.fac)
```

```
    |--------------------------------------|
    | MCM Result Summary (MCMSEM v0.25.0)  |
    |--------------------------------------|
    device         : cpu
    N phenotypes   : 5
    N latents      : 1
    Parameters summary
      label lhs edge rhs       est          se p         last_gradient
    1  a1_1  f1   =~  x1 0.7022672 0.009946645 0 -1.02152116596699e-05
    2  a1_2  f1   =~  x2 0.7064014 0.011915073 0 -1.85966491699219e-05
    3  a1_3  f1   =~  x3 0.7379681 0.010522911 0 -3.27005982398987e-05
    4  a1_4  f1   =~  x4 0.7017218 0.010998424 0   5.0276517868042e-05
    5  a1_5  f1   =~  x5 0.6471799 0.011246964 0  1.90660357475281e-05
    Variances summary
      label lhs edge rhs       est         se p         last_gradient
    1    s1  x1   ~~  x1 0.5032219 0.01076200 0   5.6014396250248e-05
    2    s2  x2   ~~  x2 0.5905972 0.01171703 0 -3.84002923965454e-05
    3    s3  x3   ~~  x3 0.4630063 0.01173199 0 -2.92276963591576e-05
    4    s4  x4   ~~  x4 0.5564786 0.01018167 0 -2.13757157325745e-05
    5    s5  x5   ~~  x5 0.6210640 0.01126981 0 -2.62558460235596e-05
    Skewness summary
      label edge v1 v2 v3        est         se            p last_gradient
    1   sk1  ~~~ x1 x1 x1 0.51329112 0.05142056 1.823493e-23             0
    2   sk2  ~~~ x2 x2 x2 0.03569529 0.05474808 5.144065e-01             0
    3   sk3  ~~~ x3 x3 x3 0.64177245 0.05309632 1.238024e-33             0
    4   sk4  ~~~ x4 x4 x4 0.27729309 0.04209400 4.473736e-11             0
    5   sk5  ~~~ x5 x5 x5 0.24263433 0.04182465 6.583126e-09             0
```
Very similar results if you compare the MCMSEM estimates to the lavaan standardized results (last column in lavaan). This inst too unexpected we fitted very similar models actually (a single factor model)
However, in MCMSEM we can actually use the multivariate skewness and kurtosis between the variables to just estimate the directed network with all paths in all directions! Let's do that now:

```
mod.network <- MCMmodel(simmdatasumm, n_latent=0,
                        causal_observed = T ,constrained_a = FALSE)
res.network <- MCMfit(mod.network , simmdatasumm,
                      optimizers=c("rprop", "lbfgs"), optim_iters=c(2000, 25),
                      learning_rate=c(0.15, 1.0), monitor_grads = TRUE)
summary(res.network)
```

```
    |--------------------------------------|
    | MCM Result Summary (MCMSEM v0.23.0)  |
    |--------------------------------------|
    device         : cpu
    N phenotypes   : 5
    N latents      : 0
    Parameters summary
       label lhs edge rhs          est         se            p         last_gradient
    1   b1_2  x1   ~>  x2  0.034553736 0.02912424 2.354544e-01 -0.000110398046672344
    2   b1_3  x1   ~>  x3  0.269620180 0.01519669 1.986098e-70  -1.9522849470377e-05
    3   b1_4  x1   ~>  x4  0.017433381 0.03216765 5.878503e-01 -0.000103633850812912
    4   b1_5  x1   ~>  x5  0.024599750 0.03638230 4.989479e-01 -0.000334467738866806
    5   b2_1  x2   ~>  x1  0.239437073 0.02666846 2.750187e-19 -6.62673264741898e-05
    6   b2_3  x2   ~>  x3  0.311603338 0.02057957 8.636872e-52 -2.61366367340088e-05
    7   b2_4  x2   ~>  x4 -0.003027891 0.05276040 9.542350e-01 -0.000234650447964668
    8   b2_5  x2   ~>  x5  0.103679739 0.05321574 5.138017e-02 -0.000123688019812107
    9   b3_1  x3   ~>  x1 -0.020487245 0.02072841 3.229738e-01 -9.44137573242188e-05
    10  b3_2  x3   ~>  x2 -0.095153973 0.02617117 2.770851e-04  5.33880665898323e-07
    11  b3_4  x3   ~>  x4 -0.027613048 0.02671749 3.013609e-01 -7.98152759671211e-05
    12  b3_5  x3   ~>  x5  0.230087474 0.02561380 2.636062e-19 -0.000242706504650414
```

Networks are really better inspected trough visualization then trough staring at path estimates, so lets go ahead and to that:

```
mod.network <- MCMmodel(simmdatasumm, n_latent=0,
                        causal_observed = T ,constrained_a = FALSE)
res.network <- MCMfit(mod.network , simmdatasumm,
                      optimizers=c("rprop", "lbfgs"), optim_iters=c(2000, 25),    
                      learning_rate=c(0.15, 1.0), monitor_grads = TRUE)

# and plot:
layout(matrix(c(1,2),1,2))
# Model
plot(res.network,layout="circle")

# Simualted Truth:
qgraph::qgraph(t(b),layout="circle",diag=T,curveAll=T)
```
Left we have the estimated network, right the true network, note that in some cases graph changes the arc of the edge but if you look at the direction you'll see these are very similar!

![Figure 1](https://raw.githubusercontent.com/zenabtamimy/MCMSEM/dev-torch/imgs/0.1.Figure1.png)

Finally, MCMSEM allows us to compare the two models in terms of fit, and in terms of fit to holdout data we generated previously.

 ```
 MCMcompareloss(list(res.single.fac,res.network),test_data = simdata.holdout)
 ```

 ```
        mse_train_loss train_chisq train_bic mse_test_loss  mse_diff mse_test_chisq mse_test_bic N_parameters
 model1    0.344228446   8605.7112 8808.2438     0.8458791        NA       8458.791     8642.998           20
 model2    0.005098701    127.4675  481.8996     0.5528885 0.2929906       5528.885     5851.246           35
 ```

So in the training data (data you used to fit the model) the loss of the network model is way lower than that of the factor model, so are the chi-square statistics, the BIC. 
in the test data we still have a lower loss for the network model, and a lower chi-square and BIC as well. The network model does have more parameters (complexity): 20 directed edges, 5 skewness parameters, 5 variances, and 5 kurtosis parameters. 
The added complexity outweighs the cost because the model does (way) better in new data. This was the advertorial, there are practical theoretical and methodological nuances and limitations, but I bet you are motivated to learn about these now!

### More information

For a more detailed description of the various functions used, see our wiki pages.

## Contribute

If you would like to contribute to MCMSEM, please do so via the [dev-torch branch](https://github.com/zenabtamimy/MCMSEM/tree/dev-torch). 
