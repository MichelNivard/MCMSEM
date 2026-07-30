# MCMSEM
R-package which allows users to run multi co-moment structural equation models.

## MCMSEM version 0.27.0
Welcome to the new and improved MCMSEM. If you want to use the MCMSEM version as it was used in [the original publication](https://doi.org/10.31235/osf.io/ynam2), please go to the [v0.1.1 release](https://github.com/zenabtamimy/MCMSEM/releases/tag/v0.1.1).

This version is considerably more powerful than our previous version. Some highlights:
 - Expanded to allow for any N variables instead of just two (use at your own risk)
 - Far more flexibility for custom model creation via MCMSEM model objects
 - More detailed fit statistics through custom MCMSEM result objects
 - Significantly improved performance, and enabled optimization on GPU
 - Asymptotic calculation of standard errors (bootstrapping no longer required)
 - Exportable data, making it easier for researchers to share moment matrices for MCMSEM without sharing raw data
 - A stationary dynamic VAR(1) kernel alongside the original contemporaneous structural kernel

If you are new to this version of MCMSEM we highly recommend reading our Wiki before starting, as the syntax for using MCMSEM has changed significantly since `v0.1.1`.

## Choosing a moment kernel

MCMSEM now makes the scientific assumptions about whether to consider, or include, time in the model as a variable explicit:

```r
# Existing behavior and the default
contemporaneous_model <- MCMmodel(ds, kernel = "contemporaneous")

# Stationary observed-state VAR(1)
dynamic_model <- MCMmodel(ds, kernel = "dynamic")
```

`kernel = "static"` is a supported, silent alias for `"contemporaneous"`.
Legacy summary and model objects without kernel metadata are also treated as
contemporaneous. The two kernels let the user specify slighlty different types of models and help answer different slightly questions.

Why would we care about this kind of nuance? We (the developers) envision people will use MCMSEM alongside other estimation methods, as the assumptions about moments are novel to many users and people would want some external validaiton. One of the methods we could see people using  MCMSEM alongside in psychology would be (random intercept) cross lagged panel models, across 3+ waves of data. Those models define "causal" paths from variabe y on x across time, where y at t-1 influences x at t, as in MCMSEM you'd model the distubances as non-guasian, we'd need a way to model the disturbance of y at t, which hasn't influences x contemoreneously at t yet. This means that to get the same estimate/estimnd out of MCMSEM as you'd get out of a stationary (RI)CLPM, you have to model explict disturbances at each time point t. This requires the new "dynamic" kernel. 

Which kernel you pick, essentially the choice you make with respect to how to model the data, depends on whether u conceive as the causal process as a static process that has unfoldeded over the past. For example a static process could be conceived of as follows: 

> taller adults are heavier, because hte volum of their bodies, all else being equal, is higher, and at similar density this means they are heavier, the process of growth is completed in adults, so there is no dynamic change in height, that results in changes in weight.

An example of a dynamic process can be found in markets (econ) or emotions (psychology), for example:

> if I dont sleep well tonight (day = t-2), il be tired tomorrow (day = t), if I then do sleep well (day = t), ill not be tired the day after tomorrow (day = t).

MCMSEM an now model cross sectional data as if its a part of a dynamic system (under assumtions like stationaiity, and non-gaussian disturbances, and all confounders being gausian.).

### Contemporaneous Structural MCMSEM

Contemporaneous structural MCMSEM assumes that the measured variables can be represented by a set of structural equations at one conceptual occasion, such as $Y=\beta X+\varepsilon_Y$. Choosing this model means treating the causal relation as meaningful without explicitly modelling the time over which it unfolds. The coefficient $\beta$ is interpreted through an intervention: replacing the equation for $X$ by $X=x$ changes the value generated for $Y$. Any prior history, adaptation, feedback, or equilibrium process is absorbed into the variables and disturbances rather than represented explicitly. Higher-order-moment identification also requires strong disturbance assumptions: the relevant structural errors must be sufficiently non-Gaussian, their dependence structure must be correctly specified, and omitted common causes must either be absent or explicitly modelled. This model is most defensible when the variables are naturally contemporaneous constructs, when one variable plausibly acts effectively before the other within the measurement window, or when the coefficient is understood as an equilibrium or total same-occasion response. It matters greatly if reciprocal processes operate within that window: a static directional path may then summarize an integrated equilibrium relationship rather than a single mechanistic transition.


The contemporaneous kernel uses the existing **Reticular Action Model (RAM)** specification. Here, $F$ maps the complete set of observed and latent variables onto the observed variables, $I$ is the identity matrix, $A$ contains the directed structural paths, and $S_2$, $S_3$, and $S_4$ contain the freely specified second-, third-, and fourth-order disturbance co-moments. The expected co-moment matrices are

$$
M_2 = F(I-A)^{-1}
S_2
(I-A)^{-T}F^\top,
$$

$$
M_3 = F(I-A)^{-1}
S_3
\left[
(I-A)^{-T}\otimes(I-A)^{-T}
\right]
\left(
F^\top\otimes F^\top
\right),
$$

and

$$
M_4 = F(I-A)^{-1}
S_4
\left[
(I-A)^{-T}\otimes
(I-A)^{-T}\otimes
(I-A)^{-T}
\right]
\left(
F^\top\otimes F^\top\otimes F^\top
\right).
$$

The factor $(I-A)^{-1}$ is the RAM total-effects matrix: it propagates each disturbance through the contemporaneous directed paths in $A$. A covariance has two indices, co-skewness has three, and co-kurtosis has four, so the same RAM transformation must act once on each index. In the matrix representation used by MCMSEM, the left-hand factor transforms the first index, while the Kronecker-product factors transform the remaining indices. The paths in $A$ are therefore estimated from the way one contemporaneous structural system jointly transforms the disturbance covariance, co-skewness, and co-kurtosis matrices. This is the notation and parameterisation already used by the original MCMSEM model.


### Stationary Dynamic MCMSEM

Stationary dynamic MCMSEM assumes instead that the world evolves through repeated, time-homogeneous transitions, $z_t=Bz_{t-1}+\varepsilon_t$. Selecting it commits you to a particular temporal resolution: the effects in $B$ occur over one chosen lag, the same transition matrix operates at every wave, the innovation distribution is stable over time, and the process has reached stationarity. It also assumes that innovations are independent across time and, for the identifying higher-order-moment argument, that the non-Gaussian innovation components have the specified independence structure. Stable Gaussian confounding may be represented separately through a residual covariance matrix. The cross-sectional distribution is then interpreted as the accumulated result of infinitely many past shocks, not as a timeless structural relation. This matters because changing the measurement interval changes the meaning and usually the numerical value of $B$: a one-day cross-lag is not the same parameter as a one-year cross-lag. The dynamic model is therefore mechanistically clearer, but it buys that clarity by imposing stronger assumptions about stationarity, lag structure, and the absence of unmodelled intermediate dynamics.


In the dynamic kernel, the matrix $B$ holds the causal paths, and it does not transform the innovations only once. Each past innovation has passed through the transition matrix a different number of times:

$$
z_t = \varepsilon_t+B\varepsilon_{t-1}+B^2\varepsilon_{t-2}+\cdots.
$$

Consequently, the stationary cumulant of order $r$ is the accumulated contribution of innovations from all previous times:

$$
(I-B^{\otimes r})^{-1}D_r.
$$

In particular,

$$
C_2=(I-B^{\otimes2})^{-1}D_2,
\qquad
C_3=(I-B^{\otimes3})^{-1}D_3,
\qquad
K_4=(I-B^{\otimes4})^{-1}D_4.
$$

The parameters in $B$ are therefore identified by the pattern produced when independent non-Gaussian innovations repeatedly propagate through the system. If correlated Gaussian residual or random-intercept components are included, their covariance $\Psi_G$ is added to $C_2$, while $C_3$ and $K_4$ remain unchanged; raw $M_4$ is then reconstructed using the total covariance $C_2+\Psi_G$.

Alternatively, `residual_family = "common_gamma"` adds one non-Gaussian
confounder to the stationary marginal state. Write
$U=(G-\alpha)/\sqrt{\alpha}$ for $G\sim\operatorname{Gamma}(\alpha,1)$ and
$y_t=z_t+\lambda U$. Then

$$
C_{2,U}=\lambda\lambda^\top,\qquad
C_{3,U}=\frac{2}{\sqrt{\alpha}}\lambda^{\otimes3},\qquad
K_{4,U}=\frac{6}{\alpha}\lambda^{\otimes4}.
$$

The loadings may be signed and `shape_Gamma` is positive. This factor replaces,
rather than supplements, the Gaussian residual in the current API. It is an
additive marginal confounder: it is not an innovation repeatedly propagated
through $B$.

```r
gamma_model <- MCMmodel(
  ds, kernel = "dynamic", residual_family = "common_gamma"
)
gamma_parameters <- MCMparameters(gamma_model)
gamma_parameters[grepl("Gamma", gamma_parameters$name), ]
```


## Free, fixed, and derived parameters

`MCMparameter()` adds auxiliary parameters or converts an existing matrix label. This can help users define parameteric constraints consitent with specific distributions (a gamma distribution for example). Only free parameters are optimized; fixed and derived parameters do not consume degrees of freedom. This signed-gamma constraint uses one positive shape per innovation instead of separate skewness and excess-kurtosis parameters:

```r
model <- MCMmodel(ds, kernel = "dynamic")
model <- MCMparameter(
  model, "shape_Earnings", "free", start = 4, transform = "positive"
)
model <- MCMparameter(
  model, "sign_Earnings", "fixed", value = -1
)
model <- MCMparameter(
  model, "tau_Earnings", "derived",
  expression = ~ sign_Earnings * 2 / sqrt(shape_Earnings)
)
model <- MCMparameter(
  model, "kappa_Earnings", "derived",
  expression = ~ 6 / shape_Earnings
)

MCMparameters(model)
MCMdegreesoffreedom(model)
```

Expressions use a small validated language rather than arbitrary R evaluation. Supported operations are `+`, `-`, `*`, `/`, `^`, `sqrt()`, `exp()`, `log()`, `softplus()`, and `logistic()`. Free starts and fitted estimates are reported on their natural scale; positive and bounded parameters use unconstrained internal optimizer coordinates. Summaries report each parameter's type and expression, derived SEs use the delta method, and fixed SEs are zero when covariance is available.

The same API works for contemporaneous models. There, diagonal `K` entries are raw standardized fourth moments, so the signed-gamma expression is
`~ 3 + 6 / shape`; dynamic `Kappa` entries are fourth cumulants and use
`~ 6 / shape`.


The current dynamic release supports observed states with at least two variables, a VAR(1) transition, fixed unit innovation variances, diagonal innovation third and fourth cumulants, and either an optional full Gaussian residual covariance or one common-gamma residual factor. It provides identity, diagonal, and full WLS moment weights plus asymptotic robust or efficient SEs. General user-defined residual distributions, latent measurement models, VAR(q), and combined contemporaneous-plus-lagged paths are not yet supported.


See
[Choosing between contemporaneous and dynamic kernels](wiki/2.3%20Choosing%20a%20kernel.md)
for the conceptual assumptions, identification conditions, continuous-time
connection, and limitations. 

The original MCM-SEM framework is described by
Tamimy et al. (2022), [*Multi Co-Moment Structural Equation Models: Discovering
Direction of Causality in the Presence of
Confounding*](https://doi.org/10.31235/osf.io/ynam2). 

The RAM specification for higher order moments was developed by: Boudt, K., Cornilly, D., & Verdonck, T. (2020). Nearest comoment estimation with unobserved factors. Journal of Econometrics, 217(2), 381–397. https://doi.org/10.1016/j.jeconom.2019.12.009

The cumulant-identification
framework and discrete Lyapunov formulation underlying the dynamic kernel build
on Cecilie Olesen Recke, Sarah Lumpp, Nataliia Kushnerchuk, Janike Oldekop,
Jiayi Li, Jane Ivy Coons, and Elina Robeva (2026), [*Identifiability in Graphical
Discrete Lyapunov Models*](https://arxiv.org/abs/2601.21818), arXiv preprint
arXiv:2601.21818.


## Citation

If you use this package please include the following citation:  
Tamimy, Z., van Bergen, E., van der Zee, M. D., Dolan, C. V., & Nivard, M. G. (2022, June 30). Multi Co-Moment Structural Equation Models: Discovering Direction of Causality in the Presence of Confounding. [https://doi.org/10.31235/osf.io/ynam2](https://doi.org/10.31235/osf.io/ynam2)

If you use the stationary dynamic kernel, also cite:

Recke, C. O., Lumpp, S., Kushnerchuk, N., Oldekop, J., Li, J., Coons, J. I., &
Robeva, E. (2026). *Identifiability in Graphical Discrete Lyapunov Models*.
arXiv preprint arXiv:2601.21818.
[https://arxiv.org/abs/2601.21818](https://arxiv.org/abs/2601.21818)

## Installation

Currently, this package is not listed on CRAN and should therefore be installed from GitHub directly.
```
library(devtools)
install_github("https://github.com/MichelNivard/MCMSEM")
```

See the wiki `Installing MCMSEM` for more details.


### MCMSEM on GPU

See the `Installing MCMSEM` wiki

## Usage

Below you will find a short rationale with usage examples of MCMSEM. For more detailed descriptions please visit the wiki.

### A rationale for MCMSEM: this is BIG, you should care!

Let's get going, in this very short pre-tutorial I'll convince you why you should read the entire tutorial, the paper(s) and consider MCMSEM for your projects. 
This is an advertorial, not a full review of the method with all its good and bad, that's left for the paper and the rest of this wiki. This is an example of the static or contemoraneous kernel:

```{r}
library(devtools)
install_github("https://github.com/MichelNivard/MCMSEM")
library(MCMSEM)
library(lavaan)

# MAke sure we might be able to replicate this
set.seed(789)
```

The basic premise is that by modeling higher order co-moments, not just covariance MCMSEM can do incredible things. 
It can (for example) estimate directional causal effects in the presence of other causal effects in the opposite direction based on continuous variables collected in a cross-sectional and observational setting.
Given certain assumptions about the confusers hold even causal effects in the presence of confounding. So to proof that to you, let me simulate data from a dense network with bidirectional causal relations between variables. 
In the simulation we only use direct paths between variables to induce correlations, no latent variables are present.

```
b <- matrix(c(   0,  .3,   0,  .1,  .15, 
                 0,   0,   0,  .1,  .15,
               .15,  .2,   0, .12,   .2,
               .15, .15,  .1,   0,   .3,
                .1, .15, .05,   0,    0), 5,5,byrow=T)
# Latent variables don't load on the indicators:
a <- matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 0, 0), ncol=2)

# use the MCMSEM internal simuation tool:
simmdata<- simulate_data(n=100000,a=a,b=b,shape=c(7, 0, 3, 4, 5), df=c(0, 8, 0, 10, 12),asdataframe = T)


#simulate holdout data from the exact same process!
simdata.holdout <-  simulate_data(n=50000,a=a,b=b,shape=c(7, 0, 3, 4, 5), df=c(0, 8, 0, 10, 12),asdataframe = T)
cor(simmdata)
```

```
          [,1]      [,2]      [,3]      [,4]      [,5]
[1,] 1.0000000 0.4423112 0.4210967 0.4450527 0.4296345
[2,] 0.4423112 1.0000000 0.4247575 0.4082202 0.3966041
[3,] 0.4210967 0.4247575 1.0000000 0.4683351 0.4801946
[4,] 0.4450527 0.4082202 0.4683351 1.0000000 0.4268477
[5,] 0.4296345 0.3966041 0.4801946 0.4268477 1.0000000
> 
```

So we generated network data, which gives rise to 5 correlated variables, we simulated 100,000 observations of 5 variables, and we made sure these variables have some skewness and kurtosis.
Then we computed the correlations between the variables, these correlations seem sort of consistent with the influence of a single latent variable (all 5 variables are correlated about equally). 
This is a known problem right? the data don't really identify a specific model, and might fit the wrong model rather well, actually. If we fit a latent variable model to these data in `lavaan` what happens?
```
### Lavaan naiveness:
# specify a single factor model:
model <- "F1 =~ V1 + V2 + V3 + V4 + V5" 

#Fit a single factor modle to the data:
single.factor.model <- sem(model,data = simmdata)

### Based on respectable fit indices, the single factor model has pretty good fit....
fitmeasures(single.factor.model,fit.measures = c("cfi","rmsea"))
```

```
  cfi rmsea 
0.992 0.046 
```

So that's a pretty solid fit to 100,000 data points which aren't normally distributed! I think many people would happily accept that the (wrong) model fits the data well in this case.
Let's look at the estimated parameters, I took liberty of omitting a part of the result (here and in other examples below) to improve readability of the page:

```
summary(single.factor.model,standardize=T)
```

```
lavaan 0.6-12 ended normally after 24 iterations

  Estimator                                         ML


Latent Variables:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
  F1 =~                                                                 
    V1                1.000                               0.786    0.656
    V2                0.887    0.003  272.326    0.000    0.697    0.626
    V3                1.091    0.004  291.463    0.000    0.858    0.689
    V4                1.421    0.005  284.822    0.000    1.117    0.666
    V5                1.293    0.005  282.941    0.000    1.017    0.659

Variances:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
   .V1                0.820    0.003  309.754    0.000    0.820    0.570
   .V2                0.756    0.002  321.343    0.000    0.756    0.608
   .V3                0.815    0.003  294.057    0.000    0.815    0.525
   .V4                1.569    0.005  305.367    0.000    1.569    0.557
   .V5                1.345    0.004  308.145    0.000    1.345    0.565
    F1                0.618    0.003  179.721    0.000    1.000    1.000
```

### Let's try MCMSEM!
Okay we can fit the exact same model in MCMSEM, this will take longer, but that's because we are using MCMSEM for a simple model while it's meant for way more complex models...

```
# Prepare the data:
simmdatasumm <- MCMdatasummary(simmdata)

# specify the MCMSEM single factor model:
mod.single.fac <- MCMmodel(simmdatasumm, n_latent=1,
                           causal_observed = FALSE, constrained_a = FALSE,
                           kernel = "contemporaneous")
res.single.fac  <- MCMfit(mod.single.fac , simmdatasumm,
                          optimizers=c("rprop", "lbfgs"), optim_iters=c(5000, 50),
                          learning_rate=c(0.15,.35), monitor_grads = TRUE,debug = T)
summary(res.single.fac)
```

```
|--------------------------------------|
| MCM Result Summary (MCMSEM v0.27.0)  |
|--------------------------------------|
device         : cpu
N phenotypes   : 5
N latents      : 1


 Parameters summary
  label lhs edge rhs       est          se p         last_gradient
1  a1_1  f1   =~  x1 0.6774216 0.002221994 0  -0.00030372804030776
2  a1_2  f1   =~  x2 0.6460607 0.002450077 0 -0.000246435403823853
3  a1_3  f1   =~  x3 0.7112606 0.002113517 0 -0.000320896506309509
4  a1_4  f1   =~  x4 0.6737278 0.002078837 0 -0.000303968787193298
5  a1_5  f1   =~  x5 0.6594946 0.002341010 0 -0.000288307666778564

Variances summary
  label lhs edge rhs       est          se p         last_gradient
1    s1  x1   ~~  x1 0.5429684 0.002322939 0 -0.000104825012385845
2    s2  x2   ~~  x2 0.7034848 0.002680491 0 -4.81307506561279e-06
3    s3  x3   ~~  x3 0.4870218 0.002574685 0 -9.28817316889763e-05
4    s4  x4   ~~  x4 0.5476845 0.002248093 0 -0.000122410012409091
5    s5  x5   ~~  x5 0.6263368 0.002268915 0 -4.79742884635925e-05

Skewness summary
  label edge v1 v2 v3        est          se             p last_gradient
1   sk1  ~~~ x1 x1 x1 0.52245688 0.009559851  0.000000e+00             0
2   sk2  ~~~ x2 x2 x2 0.02273338 0.011608215  5.018456e-02             0
3   sk3  ~~~ x3 x3 x3 0.69856018 0.010974598  0.000000e+00             0
4   sk4  ~~~ x4 x4 x4 0.27626094 0.008559417 1.535428e-228             0
5   sk5  ~~~ x5 x5 x5 0.30374101 0.008960798 7.464604e-252             0

Kurtosis summary
  label edge v1 v2 v3 v4      est         se p         last_gradient
1    k1 ~~~~ x1 x1 x1 x1 1.405791 0.02439721 0  -5.7220458984375e-06
2    k2 ~~~~ x2 x2 x2 x2 1.963156 0.03590991 0  4.76837158203125e-06
3    k3 ~~~~ x3 x3 x3 x3 1.769238 0.03414891 0   1.9073486328125e-06
4    k4 ~~~~ x4 x4 x4 x4 1.307484 0.01698693 0 -4.76837158203125e-06
5    k5 ~~~~ x5 x5 x5 x5 1.297236 0.01762279 0  -1.9073486328125e-06
```

Very similar results if you compare the MCMSEM estimates to the lavaan standardized results (last column in lavaan). This inst too unexpected we fitted very similar models actually (a single factor model).

However, in MCMSEM we can actually use the multivariate skewness and kurtosis between the variables to just estimate the directed network with all paths in all directions! Let's do that now:

```
mod.network <- MCMmodel(simmdata, n_latent=0,
                        causal_observed = TRUE, scale_data = TRUE,
                        constrained_a = FALSE,
                        kernel = "contemporaneous")
res.network <- MCMfit(mod.network , simmdata,
                      optimizers=c("rprop", "lbfgs"), optim_iters=c(5000, 50),
                      learning_rate=c(0.15,.35), monitor_grads = TRUE,debug = T)

summary(res.network)
```

```
|--------------------------------------|
| MCM Result Summary (MCMSEM v0.27.0)  |
|--------------------------------------|
device         : cpu
N phenotypes   : 5
N latents      : 0

Parameters summary
   label lhs edge rhs           est          se             p        last_gradient
1   b1_2  x1   ~>  x2  6.151331e-03 0.006816274  3.668191e-01 0.000175460241734982
2   b1_3  x1   ~>  x3  1.452449e-01 0.002979095  0.000000e+00 -4.4724001782015e-05
3   b1_4  x1   ~>  x4  1.020069e-01 0.005918295  1.429033e-66 0.000100970733910799
4   b1_5  x1   ~>  x5  6.854491e-02 0.005808502  3.866619e-32 0.000267743365839124
5   b2_1  x2   ~>  x1  2.755117e-01 0.006062157  0.000000e+00 0.000182127900188789
6   b2_3  x2   ~>  x3  1.834408e-01 0.004096662  0.000000e+00 6.35582255199552e-05
7   b2_4  x2   ~>  x4  1.017229e-01 0.009157432  1.143939e-28 0.000135798007249832
8   b2_5  x2   ~>  x5  1.106104e-01 0.009164575  1.533458e-33 0.000255572609603405
9   b3_1  x3   ~>  x1 -9.047752e-05 0.003688662  9.804310e-01 0.000118969241157174
10  b3_2  x3   ~>  x2 -3.053474e-03 0.005180531  5.555840e-01 0.000170918647199869
11  b3_4  x3   ~>  x4  7.760608e-02 0.004376185  2.302355e-70 0.000154128996655345
12  b3_5  x3   ~>  x5  4.493757e-02 0.004955088  1.201422e-19 0.000189122278243303
```

Networks are really better inspected trough visualization then trough staring at path estimates, so lets go ahead and to that:

```
# and plot:
layout(matrix(c(1,2),1,2))
# Model
plot(res.network,layout="circle")

# Simulated Truth:
b2 <- b + res.network$model$num_matrices$S
qgraph::qgraph(t(b2),layout="circle",diag=T,curveAll=T)

```
Left we have the estimated network, right the true network, note that in some cases graph changes the arc of the edge but if you look at the direction you'll see these are very similar!

![Figure 1](https://raw.githubusercontent.com/zenabtamimy/MCMSEM/dev-torch/imgs/0.1.Figure1.png)

Finally, MCMSEM allows us to compare the two models in terms of fit, and in terms of fit to holdout data we generated previously.

 ```
 MCMcompareloss(list(res.single.fac,res.network),test_data = simdata.holdout)
 ```

 ```
       mse_train_loss train_chisq   train_bic mse_test_loss  mse_diff mse_test_chisq mse_test_bic N_parameters
model1   0.4127585292 123827.5588 124079.7895    0.43718457        NA      65577.686    65816.054           20
model2   0.0003538882    106.1665    547.5703    0.01523786 0.4219467       2285.678     2702.822           35
 ```

So in the training data (data you used to fit the model) the loss of the network model is way lower than that of the factor model, so are the chi-square statistics, the BIC. In the test data we still have a lower loss for the network model, and a lower chi-square and BIC as well. The network model does have more parameters (complexity): 20 directed edges, 5 skewness parameters, 5 variances, and 5 kurtosis parameters. The added complexity outweighs the cost because the model does (way) better in new data. 

This was the advertorial, there are practical theoretical and methodological nuances and limitations, but I bet you are motivated to learn about these now!





### Worked dynamic examples

#### Controlled simulation: CLPM and dynamic MCMSEM recover the same process

First simulate 20,000 independent subjects from a stationary bivariate VAR(1).
The rows of `B_true` are current outcomes and its columns are lagged predictors.
The two innovations are mutually independent, non-Gaussian, mean zero, and
variance one—the assumptions used by the fitted dynamic kernel.

```r
library(MCMSEM)
library(lavaan)
set.seed(20260730)
n <- 20000L
B_true <- matrix(c(0.55, 0.16,
                  -0.12, 0.45), 2, 2, byrow = TRUE)
state <- matrix(0, n, 2)
panel <- array(NA_real_, dim = c(n, 2, 4))

innovation <- function(n) {
  cbind(
    rexp(n) - 1,
    (rchisq(n, df = 5) - 5) / sqrt(10)
  )
}

# Burn in for 200 transitions, then retain four consecutive waves.
for (tt in seq_len(204L)) {
  state <- state %*% t(B_true) + innovation(n)
  if (tt > 200L) panel[, , tt - 200L] <- state
}

sim_wide <- data.frame(
  X1 = panel[, 1, 1], Y1 = panel[, 2, 1],
  X2 = panel[, 1, 2], Y2 = panel[, 2, 2],
  X3 = panel[, 1, 3], Y3 = panel[, 2, 3],
  X4 = panel[, 1, 4], Y4 = panel[, 2, 4]
)
```

Fit an equality-constrained CLPM to all four waves. Because the data-generating
process has no time-invariant trait component, a standard CLPM is the matched
longitudinal estimator in this simulation.

```r
sim_clpm <- "
  X2 ~ x_ar*X1 + y_to_x*Y1
  X3 ~ x_ar*X2 + y_to_x*Y2
  X4 ~ x_ar*X3 + y_to_x*Y3
  Y2 ~ y_ar*Y1 + x_to_y*X1
  Y3 ~ y_ar*Y2 + x_to_y*X2
  Y4 ~ y_ar*Y3 + x_to_y*X3
  X1 ~~ Y1
  X2 ~~ Y2
  X3 ~~ Y3
  X4 ~~ Y4
"
sim_clpm_fit <- sem(
  sim_clpm, data = sim_wide, estimator = "MLR", meanstructure = TRUE
)
```

Now discard waves 1–3 and fit the dynamic kernel to the final marginal
cross-section only.

```r
sim_final <- sim_wide[c("X4", "Y4")]
names(sim_final) <- c("X", "Y")
sim_ds <- MCMdatasummary(
  sim_final, scale_data = FALSE, prep_asymptotic_se = TRUE,
  use_skewness = TRUE, use_kurtosis = TRUE
)
sim_model <- MCMmodel(
  sim_ds, n_latent = 0, kernel = "dynamic",
  gaussian_residual = FALSE
)
sim_model <- MCMedit(sim_model, "B", c(1, 1), "x_ar")
sim_model <- MCMedit(sim_model, "B", c(1, 2), "y_to_x")
sim_model <- MCMedit(sim_model, "B", c(2, 1), "x_to_y")
sim_model <- MCMedit(sim_model, "B", c(2, 2), "y_ar")
sim_fit <- MCMfit(
  sim_model, sim_ds, compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"), optim_iters = c(750, 40),
  learning_rate = c(0.01, 0.2),
  moment_weighting = "diagonal", se_correction = "robust",
  n_starts = 20, seed = 20260730, verbose = FALSE
)
```

The validated run recovered every transition closely:

| Path (`current <- lagged`) | Truth | CLPM estimate (SE) | Dynamic MCMSEM estimate (robust SE) |
|---|---:|---:|---:|
| X <- X | 0.550 | 0.550 (0.003) | 0.555 (0.014) |
| X <- Y | 0.160 | 0.162 (0.004) | 0.147 (0.046) |
| Y <- X | -0.120 | -0.117 (0.003) | -0.111 (0.028) |
| Y <- Y | 0.450 | 0.446 (0.004) | 0.462 (0.014) |

The CLPM had robust CFI = 1.000, RMSEA = 0.000, and SRMR = 0.003. The
dynamic solution had loss 0.000459, spectral radius 0.522, nominal df = 4,
Jacobian rank 8/8, information condition number 1.81e4, and 20/20 admissible
starts. This is the clean agreement expected when both estimators' assumptions
match the data-generating process; it is not evidence that agreement is
guaranteed with observational data.

#### Real-data illustration: earnings and hours in the 2014 SIPP panel

This example uses the U.S. Census Bureau's public-use
[2014 Survey of Income and Program Participation (SIPP) panel](https://www.census.gov/programs-surveys/sipp/data/datasets/2014-panel.html).
Its four waves cover the 2013--2016 reference years. The original Census files
are hundreds of megabytes and are not stored in this repository. Instead, the
package includes a 437 KB
[derived analysis matrix](inst/extdata/sipp_2014_panel.csv.gz) and its
[full provenance](inst/extdata/README.md).

The matrix contains only four waves of transformed log earnings and usual
weekly hours: eight columns, no identifiers, demographics, survey weights, or
source variables. The derivation selected December in each wave, retained a
working-age baseline cohort, and treated non-positive earnings or hours as
unavailable. The analysis is consequently an unweighted methodological
illustration among people with positive earnings and hours, not a
population-representative labor estimate. Variable definitions and the
public-use files are available from the Census
[wave pages](https://www.census.gov/programs-surveys/sipp/data/datasets/2014-panel/wave-1.html)
and [data dictionaries](https://www.census.gov/programs-surveys/sipp/tech-documentation/data-dictionaries/data-dictionaries-2014.html).

Load the exact matrix used below:

```r
library(MCMSEM)
library(lavaan)

sipp_file <- system.file(
  "extdata", "sipp_2014_panel.csv.gz", package = "MCMSEM"
)
if (!nzchar(sipp_file)) {
  # When running directly from a source checkout.
  sipp_file <- "inst/extdata/sipp_2014_panel.csv.gz"
}
sipp <- read.csv(sipp_file, na.strings = c("", "NA"))
stopifnot(
  identical(
    names(sipp),
    c(
      "Earnings1", "Hours1", "Earnings2", "Hours2",
      "Earnings3", "Hours3", "Earnings4", "Hours4"
    )
  )
)
```

Log earnings and hours were each centered using wave-1 complete cases and
rescaled as `2 * (value - wave1_mean) / wave1_sd`. The same affine
transformation was used at every wave and in every model. It gives both
variables variance four in wave 1, comfortably above the dynamic kernel's fixed
unit innovation variance. Thus the reported unstandardized paths share common
measurement units and are approximately standardized; scale differences do not
explain differences between methods.

Dynamic MCMSEM assumes a stationary marginal distribution, so in the
population any wave could be used. We use wave 1 because it has the largest
jointly observed sample. The empirical moments are similar but not identical,
so stationarity remains an approximation rather than something established by
choosing the largest wave:

| Wave | Reference year | Complete N | Earnings mean | Hours mean | Earnings variance | Hours variance | Covariance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2013 | 22,049 | 0.000 | 0.000 | 4.000 | 4.000 | 2.031 |
| 2 | 2014 | 15,787 | 0.036 | -0.034 | 4.192 | 4.170 | 2.179 |
| 3 | 2015 | 12,190 | 0.184 | 0.055 | 4.012 | 4.205 | 2.087 |
| 4 | 2016 | 10,446 | 0.268 | 0.033 | 4.057 | 4.156 | 2.166 |

##### Observed-score dynamics: CLPM versus plain dynamic MCMSEM

First fit a traditional CLPM. The four paths are constrained equal across the
three annual transitions. MLR and FIML retain incomplete panel rows and provide
non-normality-robust inference.

```r
clpm_syntax <- "
  Earnings2 ~ earnings_ar*Earnings1 + hours_to_earnings*Hours1
  Earnings3 ~ earnings_ar*Earnings2 + hours_to_earnings*Hours2
  Earnings4 ~ earnings_ar*Earnings3 + hours_to_earnings*Hours3

  Hours2 ~ hours_ar*Hours1 + earnings_to_hours*Earnings1
  Hours3 ~ hours_ar*Hours2 + earnings_to_hours*Earnings2
  Hours4 ~ hours_ar*Hours3 + earnings_to_hours*Earnings3

  Earnings1 ~~ Hours1
  Earnings2 ~~ Hours2
  Earnings3 ~~ Hours3
  Earnings4 ~~ Hours4
"
clpm_fit <- sem(
  clpm_syntax, data = sipp,
  estimator = "MLR", missing = "fiml", meanstructure = TRUE
)
```

Now fit the basic dynamic kernel with freely estimated innovation skewness and
kurtosis and no residual confounder. Only the largest cross-section enters
MCMSEM.

```r
complete_n <- vapply(seq_len(4L), function(wave) {
  sum(complete.cases(sipp[paste0(c("Earnings", "Hours"), wave)]))
}, integer(1))
mcm_wave <- which.max(complete_n)
mcm_values <- na.omit(
  sipp[paste0(c("Earnings", "Hours"), mcm_wave)]
)
names(mcm_values) <- c("Earnings", "Hours")

dynamic_data <- MCMdatasummary(
  mcm_values,
  scale_data = FALSE,
  prep_asymptotic_se = TRUE,
  use_skewness = TRUE,
  use_kurtosis = TRUE
)

make_dynamic_model <- function(residual_family) {
  model <- MCMmodel(
    dynamic_data,
    n_latent = 0,
    kernel = "dynamic",
    residual_family = residual_family
  )
  model <- MCMedit(model, "B", c(1, 1), "Earnings_AR")
  model <- MCMedit(
    model, "B", c(1, 2), "Hours_lag_to_Earnings"
  )
  model <- MCMedit(
    model, "B", c(2, 1), "Earnings_lag_to_Hours"
  )
  MCMedit(model, "B", c(2, 2), "Hours_AR")
}

plain_model <- make_dynamic_model("none")
plain_fit <- MCMfit(
  plain_model,
  dynamic_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(1200, 80),
  learning_rate = c(0.01, 0.005),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 40,
  seed = 20260901,
  verbose = FALSE
)
```

The validated estimates were:

| Path (current <- lagged) | CLPM estimate (robust SE) | Plain dynamic MCMSEM estimate (robust SE) |
|---|---:|---:|
| Earnings <- earnings | 0.628 (0.011) | 0.980 (0.026) |
| Earnings <- hours | 0.063 (0.008) | -0.315 (0.121) |
| Hours <- earnings | 0.124 (0.007) | 0.521 (0.112) |
| Hours <- hours | 0.476 (0.009) | 0.486 (0.113) |

The CLPM used all 24,505 contributing panel rows. Its scaled CFI, TLI, and
RMSEA were 0.905, 0.867, and 0.052; SRMR was 0.086. The plain MCMSEM solution
had loss 0.3659, spectral radius 0.800, nominal df = 4, Jacobian rank
8/8, information condition number $6.94 \times 10^6$, and 10/40 admissible
starts. The earnings autoregression is also effectively at its 0.98 upper
bound. The solution is reportable, but its boundary proximity and multistart
attrition make it a sensitivity result rather than a clean point estimate.

These are comparable observed-scale transition matrices, but they are not
estimated from the same information. The CLPM observes the transitions.
MCMSEM reconstructs a stationary transition from one marginal distribution
using non-Gaussian higher moments. The controlled simulation above shows that
they can converge when both models' assumptions hold; agreement is not
guaranteed in observational data.

##### Within-person dynamics: RI-CLPM and a Gaussian MCMSEM confounder

A CLPM mixes stable between-person differences with within-person change. The
RI-CLPM below separates two correlated random intercepts from within-person
deviations and estimates the transition matrix among those deviations.

```r
riclpm_syntax <- "
  RI_Earnings =~ 1*Earnings1 + 1*Earnings2 + 1*Earnings3 + 1*Earnings4
  RI_Hours =~ 1*Hours1 + 1*Hours2 + 1*Hours3 + 1*Hours4

  wE1 =~ 1*Earnings1
  wE2 =~ 1*Earnings2
  wE3 =~ 1*Earnings3
  wE4 =~ 1*Earnings4
  wH1 =~ 1*Hours1
  wH2 =~ 1*Hours2
  wH3 =~ 1*Hours3
  wH4 =~ 1*Hours4

  Earnings1 ~~ 0*Earnings1
  Earnings2 ~~ 0*Earnings2
  Earnings3 ~~ 0*Earnings3
  Earnings4 ~~ 0*Earnings4
  Hours1 ~~ 0*Hours1
  Hours2 ~~ 0*Hours2
  Hours3 ~~ 0*Hours3
  Hours4 ~~ 0*Hours4

  wE2 ~ earnings_ar*wE1 + hours_to_earnings*wH1
  wE3 ~ earnings_ar*wE2 + hours_to_earnings*wH2
  wE4 ~ earnings_ar*wE3 + hours_to_earnings*wH3
  wH2 ~ hours_ar*wH1 + earnings_to_hours*wE1
  wH3 ~ hours_ar*wH2 + earnings_to_hours*wE2
  wH4 ~ hours_ar*wH3 + earnings_to_hours*wE3

  RI_Earnings ~~ RI_Hours
  wE1 ~~ wH1
  wE2 ~~ wH2
  wE3 ~~ wH3
  wE4 ~~ wH4
  RI_Earnings ~~ 0*wE1 + 0*wH1
  RI_Hours ~~ 0*wE1 + 0*wH1
"
riclpm_fit <- sem(
  riclpm_syntax, data = sipp,
  estimator = "MLR", missing = "fiml",
  meanstructure = TRUE, fixed.x = FALSE
)
stopifnot(lavInspect(riclpm_fit, "post.check"))
```

In dynamic MCMSEM, `residual_family = "gaussian"` adds an unrestricted
Gaussian covariance that is not propagated through the transition matrix. It
is conceptually analogous to a joint distribution of stable between-person
differences: it can absorb variance and covariance that should not be assigned
to the dynamic innovations. Unlike the RI-CLPM, however, a one-wave MCMSEM fit
does not observe that the component persists over time; stability is a
distributional interpretation.

```r
gaussian_model <- make_dynamic_model("gaussian")
gaussian_fit <- MCMfit(
  gaussian_model,
  dynamic_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(1400, 80),
  learning_rate = c(0.01, 0.005),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 20,
  seed = 20260731,
  verbose = FALSE
)
gaussian_fit$Psi_G
```

| Path (current <- lagged) | RI-CLPM within-person estimate (robust SE) | Dynamic MCMSEM + Gaussian residual (robust SE) |
|---|---:|---:|
| Earnings <- earnings | 0.105 (0.024) | 0.627 (0.091) |
| Earnings <- hours | 0.051 (0.012) | 0.127 (0.104) |
| Hours <- earnings | 0.009 (0.015) | 0.465 (0.116) |
| Hours <- hours | 0.165 (0.017) | 0.591 (0.163) |

The RI-CLPM passed lavaan's post-estimation check and fit the longitudinal
covariance structure closely: scaled CFI = 0.998, TLI = 0.996, RMSEA = 0.009,
and SRMR = 0.012. The MCMSEM Gaussian residual covariance was

```text
          Earnings  Hours
Earnings     1.842  0.450
Hours        0.450  0.421
```

which implies a residual correlation of 0.511. Its loss was 0.1004, spectral
radius 0.852, nominal df = 1, and Jacobian rank 11/11. The information
condition number was $1.17 \times 10^8$, and the robust SEs above use the
ordinary un-ridged sandwich information matrix. The full-rank Jacobian does not
make this nearly saturated decomposition precise; the large condition number
and wide SEs are central results.

##### Could the confounder itself be non-Gaussian?

A Gaussian residual affects covariance and Gaussian fourth-moment pairings but
has no third or fourth cumulants. If the stable source is skewed, assigning all
higher cumulants to the dynamic innovations may distort the transition matrix.
The common-gamma family adds one centered variance-one gamma factor `U` with
signed loadings `lambda`:

$$
X_t = X_t^{\mathrm{dynamic}} + \lambda U,
\qquad
\operatorname{Cov}(\lambda U)=\lambda\lambda^\prime.
$$

Its positive shape `alpha` determines skewness $2/\sqrt{\alpha}$ and
excess kurtosis $6/\alpha$. Large `alpha` approaches a rank-one Gaussian
factor. Loadings remain signed, so odd cumulants retain their direction.

The direct API is:

```r
gamma_model <- make_dynamic_model("common_gamma")
```

With freely estimated innovation skewness and kurtosis, this bivariate model
has only one overidentifying degree of freedom and admits competing,
ill-conditioned decompositions. For the reported test we therefore imposed a
scientifically explicit signed-gamma relationship on each innovation:
earnings innovations were allowed negative skew and hours innovations positive
skew, while their two shapes remained free. For example:

```r
gamma_model <- MCMparameter(
  gamma_model, "shape_Earnings", "free",
  start = 0.114, transform = "positive"
)
gamma_model <- MCMparameter(
  gamma_model, "shape_Hours", "free",
  start = 0.295, transform = "positive"
)
gamma_model <- MCMparameter(
  gamma_model, "sign_Earnings", "fixed", value = -1
)
gamma_model <- MCMparameter(
  gamma_model, "sign_Hours", "fixed", value = 1
)
gamma_model <- MCMparameter(
  gamma_model, "tau_Earnings", "derived",
  expression = ~ sign_Earnings * 2 / sqrt(shape_Earnings)
)
gamma_model <- MCMparameter(
  gamma_model, "kappa_Earnings", "derived",
  expression = ~ 6 / shape_Earnings
)
gamma_model <- MCMparameter(
  gamma_model, "tau_Hours", "derived",
  expression = ~ sign_Hours * 2 / sqrt(shape_Hours)
)
gamma_model <- MCMparameter(
  gamma_model, "kappa_Hours", "derived",
  expression = ~ 6 / shape_Hours
)

# A direct multistart fit. The validation script below uses a broader explicit
# grid over confounder shapes and loading orientations.
gamma_fit <- MCMfit(
  gamma_model,
  dynamic_data,
  compute_se = TRUE,
  optimizers = c("rprop", "lbfgs"),
  optim_iters = c(1200, 100),
  learning_rate = c(0.01, 0.005),
  moment_weighting = "diagonal",
  se_correction = "robust",
  n_starts = 30,
  seed = 20260902,
  verbose = FALSE
)
gamma_fit$dynamic$common_gamma
```

The validation script searches both loading orientations and starting shapes
0.25, 2, and 20 before computing robust SEs. That search gave:

| Path (current <- lagged) | Common-gamma MCMSEM estimate (robust SE) |
|---|---:|
| Earnings <- earnings | 0.676 (0.032) |
| Earnings <- hours | 0.075 (0.052) |
| Hours <- earnings | 0.422 (0.030) |
| Hours <- hours | 0.651 (0.036) |

The common-factor loadings were -1.347 (SE 0.108) for earnings and -0.360
(SE 0.235) for hours. Their product gives a positive rank-one covariance. The
estimated gamma shape was 153.8 (SE 398.2), implying skewness 0.161
(SE 0.209) and excess kurtosis 0.039 (SE 0.101). The loss was 0.1090 and the
information condition number was $7.29 \times 10^{10}$.

This fit does not provide reliable evidence that the confounder is
non-Gaussian: the point estimate is close to the Gaussian limit, its
distributional SEs are large, and the overall decomposition is extremely
ill-conditioned. That conclusion is conditional on the signed-gamma innovation
constraints; it is not a general test that every possible confounder is
Gaussian.

That qualification matters empirically. An optional grid with the innovation
third and fourth cumulants left free found another basin with loss 0.00024,
shape 0.489 (SE 0.087), and loadings 1.213 (SE 0.039) and -0.064
(SE 0.037). Its information condition was $1.30 \times 10^8$ and it had only
one nominal df. This solution describes a strongly non-Gaussian component
almost entirely specific to earnings, not a convincing shared earnings-hours
confounder. Run that longer sensitivity analysis with

```sh
Rscript inst/validation/longitudinal_clpm_example.R --unrestricted-gamma
```

The constrained and unrestricted results together show that the data do not
support a specification-invariant conclusion about confounder shape.

##### How to interpret agreement and disagreement

The three comparisons should be read as triangulation rather than as competing
software implementations of one regression:

1. The CLPM and plain dynamic MCMSEM both describe observed-score dynamics.
   The former identifies paths from repeated transitions; the latter identifies
   them from a stationary marginal distribution and higher cumulants.
2. The RI-CLPM and residual-adjusted MCMSEM both attempt a within-versus-between
   decomposition. The RI-CLPM observes stable components across waves; MCMSEM
   infers a residual distribution under stronger assumptions.
3. A residual family can improve the scientific match while weakening
   identification. Compare Jacobian rank, information condition, nominal df,
   spectral radius, multistart behavior, and robust SEs—not only training loss.
4. The controlled simulation demonstrates that CLPM and MCMSEM can converge to
   the same transition matrix. Real-data estimates do not have to converge:
   disagreement can reflect different estimands, nonstationarity, weak
   distributional identification, or model misspecification rather than a
   scaling error.

The complete executable analysis, including the simulation, all multistart
settings, wave diagnostics, and machine-readable output, is in
[`inst/validation/longitudinal_clpm_example.R`](inst/validation/longitudinal_clpm_example.R);
a compact record of the validated results is in
[`inst/validation/longitudinal_clpm_example.md`](inst/validation/longitudinal_clpm_example.md).

### More information

For a more detailed description of the various functions used, see our wiki pages.

## Contribute

If you would like to contribute to MCMSEM, please do so via the [dev-torch branch](https://github.com/zenabtamimy/MCMSEM/tree/dev-torch). 
