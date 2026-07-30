# Longitudinal models and largest-wave dynamic MCMSEM

The validation script runs a controlled simulation followed by an unweighted
real-data illustration using a minimal derived SIPP analysis matrix:

```sh
Rscript inst/validation/longitudinal_clpm_example.R
```

Generated CSV summaries are written below
`validation-output/longitudinal-example/`, which is excluded from Git and
source packages.

## Controlled simulation without a stable confounder

The simulation generates 20,000 independent subjects, burns in a bivariate
VAR(1) for 200 transitions, and retains four consecutive waves. The transition
matrix is

```text
          lagged X  lagged Y
current X     0.55      0.16
current Y    -0.12      0.45
```

The first innovation is centered exponential and the second centered,
variance-one chi-square with five degrees of freedom. An
equality-constrained CLPM uses all four waves. Plain dynamic MCMSEM receives
only the final marginal cross-section and uses diagonal WLS, robust sandwich
SEs, and 20 starts.

| Path (current <- lagged) | Truth | CLPM estimate (SE) | Dynamic MCMSEM estimate (robust SE) |
|---|---:|---:|---:|
| X <- X | 0.550 | 0.550 (0.003) | 0.555 (0.014) |
| X <- Y | 0.160 | 0.162 (0.004) | 0.147 (0.046) |
| Y <- X | -0.120 | -0.117 (0.003) | -0.111 (0.028) |
| Y <- Y | 0.450 | 0.446 (0.004) | 0.462 (0.014) |

The CLPM's robust CFI, TLI, RMSEA, and SRMR were 1.000, 1.000, 0.000,
and 0.003. The dynamic fit had loss 0.000459, spectral radius 0.522,
nominal df = 4, Jacobian rank 8/8, information condition
$1.81 \times 10^4$, and 20/20 admissible starts. This is the expected
calibration result when the two estimators' assumptions and estimands match.

## Controlled simulation with a Gaussian stable confounder

A second simulation generates 100,000 subjects from the same transition
matrix. Its independent innovations are centered, variance-one gamma variables
with shapes 1 and 2.5. A time-invariant bivariate Gaussian random intercept
with covariance

```text
      X    Y
X  0.50 0.20
Y  0.20 0.35
```

is added to every wave. This simultaneously satisfies the RI-CLPM's stable
random-intercept assumptions and dynamic MCMSEM's additive Gaussian residual
assumptions. The MCMSEM fit constrains each innovation's skewness and excess
kurtosis to its freely estimated positive shape, leaving three nominal df.
One of ten starts is initialized at the known generating values; the other
nine are randomized.

| Path (current <- lagged) | Truth | CLPM estimate (SE) | RI-CLPM estimate (SE) | Gaussian MCMSEM estimate (robust SE) |
|---|---:|---:|---:|---:|
| X <- X | 0.550 | 0.646 (0.001) | 0.548 (0.003) | 0.530 (0.023) |
| X <- Y | 0.160 | 0.172 (0.001) | 0.165 (0.003) | 0.207 (0.040) |
| Y <- X | -0.120 | -0.046 (0.001) | -0.118 (0.003) | -0.101 (0.016) |
| Y <- Y | 0.450 | 0.576 (0.002) | 0.452 (0.003) | 0.476 (0.028) |

All generating paths are inside the MCMSEM robust 95% intervals. Its estimated
Gaussian residual covariance was

```text
      X     Y
X 0.485 0.132
Y 0.132 0.331
```

MCMSEM had loss 0.000222, spectral radius 0.523, rank 9/9, and information
condition $2.08 \times 10^4$. The RI-CLPM had robust CFI = 1.000,
RMSEA = 0.001, and SRMR = 0.002.

## Derived SIPP matrix

The real-data illustration uses the U.S. Census Bureau's public-use 2014
Survey of Income and Program Participation panel, waves 1--4, covering the
2013--2016 reference years. The package contains only
`inst/extdata/sipp_2014_panel.csv.gz`: a 437 KB matrix with four waves of
transformed log earnings and hours. It contains 24,505 contributing rows,
eight columns, and no identifiers, demographics, weights, or Census source
variables.

The exact official sources, extraction filters, transformation constants,
sample counts, and checksum are recorded in `inst/extdata/README.md`. The
large Census public-use source files are not included. Because survey weights
are omitted, this is a methodological illustration rather than a
population-representative labor analysis.

Both variables use the same wave-1 reference transformation at every wave:
`2 * (value - wave1_mean) / wave1_sd`. Consequently, paths in all models
share measurement units and are approximately standardized.

| Wave | Year | Complete N | Earnings mean | Hours mean | Earnings variance | Hours variance | Covariance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2013 | 22,049 | 0.000 | 0.000 | 4.000 | 4.000 | 2.031 |
| 2 | 2014 | 15,787 | 0.036 | -0.034 | 4.192 | 4.170 | 2.179 |
| 3 | 2015 | 12,190 | 0.184 | 0.055 | 4.012 | 4.205 | 2.087 |
| 4 | 2016 | 10,446 | 0.268 | 0.033 | 4.057 | 4.156 | 2.166 |

Wave 1 supplies MCMSEM because it has the largest complete bivariate sample.
The modest drift makes this a precision choice under approximate—not
established—stationarity.

## SIPP: CLPM, RI-CLPM, and Gaussian-residual dynamic MCMSEM

The RI-CLPM separates correlated stable random intercepts from within-person
deviations. Dynamic MCMSEM's full Gaussian residual covariance is conceptually
analogous to a joint stable distribution, but a one-wave fit assumes rather
than observes its temporal stability.

| Path (current <- lagged) | CLPM estimate (robust SE) | RI-CLPM estimate (robust SE) | Gaussian MCMSEM estimate (robust SE) |
|---|---:|---:|---:|
| Earnings <- earnings | 0.628 (0.011) | 0.105 (0.024) | 0.627 (0.091) |
| Earnings <- hours | 0.063 (0.008) | 0.051 (0.012) | 0.127 (0.104) |
| Hours <- earnings | 0.124 (0.007) | 0.009 (0.015) | 0.465 (0.116) |
| Hours <- hours | 0.476 (0.009) | 0.165 (0.017) | 0.591 (0.163) |

The CLPM's scaled CFI/TLI/RMSEA were 0.905/0.867/0.052 and SRMR was
0.086. The RI-CLPM passed lavaan's post-estimation check; scaled
CFI/TLI/RMSEA were 0.998/0.996/0.009 and SRMR was 0.012. The MCMSEM
Gaussian residual covariance was

```text
          Earnings  Hours
Earnings     1.842  0.450
Hours        0.450  0.421
```

and its implied correlation was 0.511. The dynamic fit had loss 0.1004,
spectral radius 0.852, df = 1, rank 11/11, and information condition
$1.17 \times 10^8$. Its full-rank but weakly conditioned decomposition and
wide un-ridged sandwich SEs should be emphasized over the small training loss.

## Common-gamma confounder

`residual_family = "common_gamma"` replaces the Gaussian residual vector
with one centered, variance-one gamma source, signed loadings, and a positive
shape. Shape determines confounder skewness and excess kurtosis. Leaving the
diagonal innovation third and fourth cumulants free gives a gamma residual or
confounder with otherwise unspecified skewed and kurtotic innovations. This is
a reasonable specification through the fitted fourth order; innovations still
have fixed unit variances, mutual independence, and diagonal higher cumulants.
The bivariate model has only one nominal df
and competing weakly identified decompositions. Constraining each innovation's
skewness and kurtosis to a signed-gamma relationship leaves three nominal df
and provides a useful more restrictive sensitivity analysis.

| Path (current <- lagged) | CLPM (robust SE) | RI-CLPM (robust SE) | Gaussian MCMSEM (robust SE) | Common gamma + signed-gamma innovations (robust SE) | Common gamma + free innovation cumulants (robust SE) |
|---|---:|---:|---:|---:|---:|
| Earnings <- earnings | 0.628 (0.011) | 0.105 (0.024) | 0.627 (0.091) | 0.676 (0.032) | 0.494 (0.058) |
| Earnings <- hours | 0.063 (0.008) | 0.051 (0.012) | 0.127 (0.104) | 0.075 (0.052) | 0.284 (0.047) |
| Hours <- earnings | 0.124 (0.007) | 0.009 (0.015) | 0.465 (0.116) | 0.422 (0.030) | 0.459 (0.046) |
| Hours <- hours | 0.476 (0.009) | 0.165 (0.017) | 0.591 (0.163) | 0.651 (0.036) | 0.579 (0.037) |

The residual loadings were -1.347 (SE 0.108) and -0.360 (SE 0.235).
The common shape was 153.8 (SE 398.2), corresponding to skewness 0.161
(SE 0.209) and excess kurtosis 0.039 (SE 0.101). Loss was 0.1090 and the
information condition was $7.29 \times 10^{10}$. Under the stated
innovation constraints, this fit does not provide reliable evidence of a
non-Gaussian confounder: the point is near the Gaussian limit and its
distributional uncertainty is large.

The longer optional unrestricted-innovation grid can be run with

```sh
Rscript inst/validation/longitudinal_clpm_example.R --unrestricted-gamma
```

It found a basin with loss 0.000231, shape 0.486 (SE 0.086), and loadings 1.214
(SE 0.039) for earnings and -0.063 (SE 0.037) for hours. Its information
condition was $1.30 \times 10^8$, its Jacobian rank was 11/11, and it had one
nominal df. The transition estimates were 0.494 for earnings autoregression,
0.579 for hours autoregression, 0.284 for hours to later earnings, and 0.459
for earnings to later hours, with robust SEs 0.058, 0.037, 0.047, and 0.046,
respectively. Thus the innovation specification matters substantively.
The free innovation third/fourth cumulants were -10.943 (SE 1.870) and 73.207
(SE 9.859) for earnings, and 4.233 (SE 0.463) and 24.045 (SE 1.887) for hours.
With fixed unit innovation variances, these are innovation skewness and excess
kurtosis, so this basin retains strongly non-Gaussian innovations as well as a
gamma residual component.
Several positive-loading-orientation starts converged near this basin, whereas
the negative-orientation fits had appreciably higher loss. This supports the
basin's reproducibility but does not resolve its one-df weak identification.
Although this model permits a shared gamma confounder, the fitted gamma
component is effectively earnings-specific because the hours loading is near
zero. It is not robust evidence for a shared earnings-hours confounder, and the
conclusion about confounder shape is constraint-dependent.

## Interpretation

- The simulations demonstrate that CLPM and plain MCMSEM can agree without a
  stable component, and that RI-CLPM and Gaussian-residual MCMSEM can agree
  when a Gaussian stable component is present.
- CLPM identifies paths from repeated transitions; MCMSEM reconstructs them
  from one stationary marginal distribution and higher cumulants.
- RI-CLPM and residual-adjusted MCMSEM both attempt a within/between
  decomposition. Only RI-CLPM directly observes persistence across waves.
- Real-data estimates need not converge. Differences may reflect estimands,
  nonstationarity, weak identification, or misspecification rather than scale.
- Training loss is insufficient. Report nominal df, Jacobian rank and
  condition, spectral radius, bounds, multistart behavior, and robust SEs.
