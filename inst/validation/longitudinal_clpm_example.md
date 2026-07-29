# Longitudinal models and largest-wave dynamic MCMSEM

The validation script runs two analyses: a controlled simulation in which both
estimators target the known VAR(1) transition matrix, followed by a real-data
sensitivity analysis.

The exact analyses can be reproduced with:

```sh
Rscript inst/validation/longitudinal_clpm_example.R
```

## Controlled stationary simulation

The simulation generates 20,000 independent subjects, burns in a bivariate
VAR(1) for 200 transitions, and retains four consecutive waves. The true
transition matrix is

```text
          lagged X  lagged Y
current X     0.55      0.16
current Y    -0.12      0.45
```

The first innovation is a centered unit-rate exponential variable and the
second is a centered and standardized chi-square variable with five degrees of
freedom. Thus both innovations have mean zero and variance one, are mutually
independent and non-Gaussian, and have distinct higher-order cumulants.

An equality-constrained CLPM uses all four waves. Dynamic MCMSEM receives only
the final marginal cross-section, uses no Gaussian residual, diagonal WLS,
robust sandwich SEs, and 20 starts.

| Path (`current <- lagged`) | Truth | CLPM estimate (SE) | Dynamic MCMSEM estimate (robust SE) |
|---|---:|---:|---:|
| X <- X | 0.550 | 0.550 (0.003) | 0.555 (0.014) |
| X <- Y | 0.160 | 0.162 (0.004) | 0.147 (0.046) |
| Y <- X | -0.120 | -0.117 (0.003) | -0.111 (0.028) |
| Y <- Y | 0.450 | 0.446 (0.004) | 0.462 (0.014) |

The CLPM's robust CFI, TLI, RMSEA, and SRMR were 1.000, 1.000, 0.000, and
0.003. The dynamic solution had loss 0.000459, spectral radius 0.522, nominal
df = 4, Jacobian rank 8/8, information condition number 1.81e4, and 20/20
admissible starts. The close agreement is the expected calibration result when
the estimators' assumptions match the data-generating process.

## Real-data sensitivity analysis

This validation uses the public [`nlswork`
extract](https://vincentarelbundock.github.io/Rdatasets/doc/sampleSelection/nlswork.html)
from the [National Longitudinal Survey of Young
Women](https://www.nlsinfo.org/content/cohorts/young-women). The data contain
28,534 observations from 4,711 women interviewed between 1968 and 1988. The
comparison uses four equally spaced two-year waves: 1971, 1973, 1975, and 1977.

The script downloads the validated [Rdatasets CSV](https://vincentarelbundock.github.io/Rdatasets/csv/sampleSelection/nlswork.csv)
(MD5 `f546ffe0bee86acb5d79b8775d341709`). It uses `4 * ln_wage` and
`hours / 5` in both analyses. These transparent linear transformations put the
observed variances on a scale compatible with the dynamic kernel's fixed unit
innovation variances; they are not z-scores.

## Models

The longitudinal analyses are a traditional CLPM and an RI-CLPM fitted with
robust maximum likelihood and FIML. Autoregressive and cross-lagged paths are
constrained equal over the three two-year transitions. The RI-CLPM separates
stable between-person intercepts from within-person deviations. Dynamic
MCMSEM is fitted only to the wave with the most complete bivariate
observations. That is 1977, with 2,167 cases. Both dynamic models use diagonal
WLS with robust sandwich SEs; the eight-parameter model has no Gaussian
residual and uses 20 starts, while the 11-parameter Gaussian-residual model
uses 30 starts.

Under exact stationarity, any wave has the same population marginal
distribution, so choosing the largest complete wave improves precision without
changing the target. The empirical wave diagnostics are not identical:

| Year | Complete N | Wage mean | Hours mean | Wage variance | Hours variance | Covariance |
|---:|---:|---:|---:|---:|---:|---:|
| 1971 | 1,851 | 6.187 | 7.331 | 2.748 | 3.613 | 0.181 |
| 1973 | 1,981 | 6.314 | 7.218 | 2.955 | 4.040 | 0.260 |
| 1975 | 2,131 | 6.327 | 7.340 | 2.650 | 3.581 | 0.010 |
| 1977 | 2,167 | 6.637 | 7.222 | 2.973 | 3.985 | 0.262 |

The larger 1977 sample is therefore used, but the drift—especially in wage
means and covariance—is evidence that stationarity is only approximate here.

## Results from the validated run

The CLPM used 3,314 participants with at least one observation. Its robust fit
indices were CFI = 0.929, TLI = 0.901, RMSEA = 0.095, and SRMR = 0.058.

| CLPM path | Estimate | Robust SE | p-value | 95% CI |
|---|---:|---:|---:|---:|
| Wage autoregression | 0.681 | 0.018 | <0.001 | [0.646, 0.715] |
| Hours to later wage | 0.039 | 0.012 | 0.001 | [0.015, 0.063] |
| Hours autoregression | 0.443 | 0.023 | <0.001 | [0.397, 0.488] |
| Wage to later hours | 0.035 | 0.020 | 0.071 | [-0.003, 0.074] |

The RI-CLPM passed `lavaan`'s post-estimation check. Its robust CFI, TLI,
RMSEA, and SRMR were 0.969, 0.949, 0.068, and 0.044.

| RI-CLPM within-person path | Estimate | Robust SE | p-value | 95% CI |
|---|---:|---:|---:|---:|
| Wage autoregression | 0.401 | 0.062 | <0.001 | [0.280, 0.522] |
| Hours to later wage | 0.056 | 0.025 | 0.026 | [0.007, 0.105] |
| Hours autoregression | 0.315 | 0.051 | <0.001 | [0.215, 0.415] |
| Wage to later hours | 0.077 | 0.060 | 0.203 | [-0.041, 0.194] |

| Dynamic MCMSEM path | Estimate | Robust SE | p-value | 95% CI |
|---|---:|---:|---:|---:|
| Wage autoregression | 0.730 | 0.088 | <0.001 | [0.558, 0.902] |
| Hours to later wage | 0.273 | 0.157 | 0.082 | [-0.035, 0.581] |
| Hours autoregression | 0.821 | 0.075 | <0.001 | [0.673, 0.969] |
| Wage to later hours | -0.391 | 0.239 | 0.101 | [-0.858, 0.077] |

The dynamic solution had spectral radius 0.840, nominal df = 4, Jacobian rank
8/8, information condition number 2.51e6, and 18/20 admissible starts. The
best diagonal-WLS loss was 0.7170.

Allowing a full Gaussian residual covariance produced:

| Dynamic MCMSEM path with Gaussian residual | Estimate | Robust SE | p-value | 95% CI |
|---|---:|---:|---:|---:|
| Wage autoregression | 0.692 | 0.402 | 0.085 | [-0.097, 1.481] |
| Hours to later wage | 0.307 | 0.124 | 0.013 | [0.064, 0.549] |
| Hours autoregression | 0.750 | 0.201 | <0.001 | [0.356, 1.145] |
| Wage to later hours | -0.518 | 0.517 | 0.316 | [-1.532, 0.495] |

Its estimated Gaussian covariance was
`matrix(c(0.343, 0.215, 0.215, 0.135), 2, 2)` and was numerically close to
rank one. Loss fell to 0.6932, but the
model had only one nominal df, information condition number `1e8`, and 20/30
admissible starts. The extra covariance decomposition is therefore weakly
identified here; the apparently lower loss should not override its much wider
SEs and poorer conditioning.

This example is a workflow comparison, not a claim that the four analyses have
the same target. The CLPM conditions on repeated individual measurements. The
dynamic MCMSEM infers a stationary transition matrix from one marginal
cross-section under independent non-Gaussian innovations and fixed innovation
scale. The different cross-lag estimates, middling CLPM fit, and high dynamic
information condition number are reasons to report diagnostics and avoid using
agreement between methods as an automatic validity criterion.

## Practical evaluation

- Decide whether the longitudinal target is observed-score dynamics (CLPM) or
  within-person deviations after stable traits are separated (RI-CLPM).
- Treat `gaussian_residual = TRUE` as a scientific model for an independent
  stable Gaussian component, not as an automatic fit improvement. In two
  variables it adds three parameters and reduces nominal df from four to one.
- Compare held-out moment loss where feasible, Jacobian rank, information
  condition, spectral radius, bounds, admissible starts, and path stability.
  Do not choose only on training loss.
- Examine marginal moments across candidate waves. Exact stationarity makes
  wave choice irrelevant in the population; empirical drift makes the
  largest-N wave a precision choice under an approximation, not proof of
  interchangeability.
