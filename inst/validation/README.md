# Opt-in methodological validation

These scripts are opt-in methodological checks and are not run during package
tests. Generated `validation-output/` files are deliberately excluded from
source builds and Git.

## Longitudinal CLPM/RI-CLPM comparison

`longitudinal_clpm_example.R` first runs both controlled stationary VAR(1)
simulations reported in the README: one without a stable component, and one
with a time-invariant bivariate Gaussian random intercept that matches both an
RI-CLPM and Gaussian-residual dynamic MCMSEM. It then reads the bundled 437 KB
`inst/extdata/sipp_2014_panel.csv.gz` analysis matrix and fits an
equality-constrained CLPM, an RI-CLPM, and dynamic MCMSEM specifications with
a full Gaussian residual covariance and a common-gamma
confounder. The SIPP analysis uses diagonal WLS with robust sandwich SEs and
explicit multistart searches.

The bundled matrix has only four waves of transformed log earnings and hours.
It contains no identifiers, demographics, survey weights, or Census source
columns. Its precise public-use source, filters, transformation, counts, and
checksum are recorded in `inst/extdata/README.md`; the large Census files are
not vendored. Run the complete opt-in analysis from the package source tree:

```sh
Rscript inst/validation/longitudinal_clpm_example.R
```

To run only the two controlled calibrations:

```sh
Rscript inst/validation/longitudinal_clpm_example.R --simulation-only
```

The main common-gamma analysis uses signed-gamma innovation constraints to
retain three nominal overidentifying df. A longer optional grid leaves the
innovation third/fourth cumulants unrestricted and documents the competing
one-df basin reported in the README:

```sh
Rscript inst/validation/longitudinal_clpm_example.R --unrestricted-gamma
```

This is an unweighted methodological illustration rather than a substantive
population estimate. The MCMSEM fits deliberately use wave 1 because it has
the largest jointly observed sample (`N = 22,049`); the script reports moments
at every wave so that the stationarity approximation remains visible.

## Common-gamma confounder validation

`common_gamma_monte_carlo.R` simulates signed-gamma innovations in a stationary
bivariate VAR(1), adds one centered common-gamma factor with signed loadings,
fits the matched constrained model, and compares empirical parameter variation
with robust delta-method SEs. It also checks the three derived entries of the
rank-one residual covariance. The default uses 50 replications:

```sh
Rscript inst/validation/common_gamma_monte_carlo.R
```

A fast two-replication smoke run is documented at the top of the script.

The following 50-replication design was run on 2026-07-29. One start was
initialized at the generating values to isolate local estimator/SE behavior
from global start selection:

```sh
Rscript inst/validation/common_gamma_monte_carlo.R \
  --repetitions=50 --n=5000 --starts=1 \
  --rprop-iters=120 --lbfgs-iters=6 \
  --output=validation-output/common-gamma-50
```

All 50 fits were stationary and had Jacobian rank 9/9. The median information
condition number was `2.33e6` (range `2.18e5` to `3.84e8`). The bivariate model
has only three nominal overidentifying df, and its robust Wald SEs were
conservative, especially for transition and factor-loading parameters. This is
a weak-identification stress test, not a demonstration of excellent Wald
calibration.

| Parameter | Bias | Empirical SD | Median ASE | Median ASE / SD | Coverage |
|---|---:|---:|---:|---:|---:|
| `phi_X` | -0.0507 | 0.0991 | 0.5716 | 5.77 | 0.98 |
| `phi_Y` | -0.0357 | 0.0924 | 0.2031 | 2.20 | 0.96 |
| `X_lag_to_Y` | 0.0142 | 0.1254 | 0.4569 | 3.64 | 0.98 |
| `Y_lag_to_X` | -0.0630 | 0.1854 | 0.3021 | 1.63 | 1.00 |
| `loading_Gamma_X` | -0.0057 | 0.1207 | 0.4223 | 3.50 | 1.00 |
| `loading_Gamma_Y` | 0.0252 | 0.0986 | 0.2048 | 2.08 | 0.94 |
| `shape_Gamma` | 0.6185 | 3.6754 | 4.7991 | 1.31 | 0.92 |
| `tau_X` | -0.0474 | 0.1711 | 0.1871 | 1.09 | 0.98 |
| `tau_Y` | 0.0121 | 0.0704 | 0.0892 | 1.27 | 0.98 |

The exact population, base/Torch, gradient, permutation, and delta-method tests
remain the primary correctness checks. For empirical work, increase
overidentification where scientifically defensible and inspect the information
condition and multistart distribution; full Jacobian rank alone is not enough.

## Standard-error calibration recorded on 2026-07-29

The command below simulated 50 stationary cross-sections with `N = 5000`, used
diagonal WLS and the robust sandwich covariance, and initialized optimization at
the data-generating values. This isolates first-order SE calibration from
multi-start search failures.

```sh
Rscript inst/validation/dynamic_se_monte_carlo.R \
  --repetitions=50 --n=5000 --starts=1 \
  --rprop-iters=150 --lbfgs-iters=8 \
  --moment-weighting=diagonal --se-correction=auto \
  --gaussian-residual=false
```

All 50 fits of the eight-parameter model without an additive Gaussian residual
were stationary, admissible, and had Jacobian rank 8. Mean ASE divided by
empirical Monte Carlo SD was close to one for both autoregressions, both
cross-lags, and both third cumulants. Fourth-cumulant calibration was less
precise.

| Parameter | Empirical SD | Mean ASE | ASE / SD | 95% coverage |
|---|---:|---:|---:|---:|
| `phi_X` | 0.0698 | 0.0746 | 1.069 | 0.98 |
| `phi_Y` | 0.0701 | 0.0744 | 1.061 | 0.98 |
| `X_lag_to_Y` | 0.1060 | 0.1033 | 0.975 | 0.90 |
| `Y_lag_to_X` | 0.1297 | 0.1371 | 1.057 | 0.92 |
| `tau_X` | 0.0983 | 0.1069 | 1.087 | 0.98 |
| `tau_Y` | 0.1012 | 0.1033 | 1.021 | 0.98 |
| `kappa_X` | 0.4771 | 0.3626 | 0.760 | 0.92 |
| `kappa_Y` | 0.6108 | 0.5895 | 0.965 | 0.86 |

The same 50-replication design with a free three-parameter Gaussian residual
covariance is nearly saturated: 12 moments and 11 parameters. It produced 47
admissible fits, a typical smallest unscaled Jacobian singular value near
`0.006`, and a weighted information condition number around `3.7e5` in a
representative sample. Unconstrained sandwich ASEs were strongly conservative
relative to the constrained Monte Carlo distribution (ASE/SD ratios from about
1.4 to 20). This is retained as a weak-identification stress test, not evidence
of acceptable calibration. Users should inspect Jacobian rank and conditioning,
use adequate overidentification, and treat Wald inference cautiously when
bounds or stationarity constraints are active.

The default script retains the Gaussian residual to exercise delta-method
`Psi_G` SEs. Pass `--gaussian-residual=false` for the better-conditioned
calibration design above.

## Free/fixed/derived parameter validation

`parameter_constraints_validation.R` is the opt-in validation for the general
parameter graph and signed-gamma demonstration in both kernels. It uses only
synthetic data. The default run uses 50 replications per kernel, a fixed master
seed, constrained and unconstrained fits with identical generating values, a
deliberately wrong common-shape constraint, and non-gamma mixture innovations:

```sh
Rscript inst/validation/parameter_constraints_validation.R \
  --repetitions=50 --n=2000 --master_seed=20260729 \
  --rprop_iters=80 --lbfgs_iters=5 --starts=2 \
  --output=validation-output/parameter-constraints
```

It writes compact CSV files for exact cumulant identities, recovery, Monte
Carlo SE/coverage calibration (including Wilson uncertainty intervals),
constraint comparisons, and misspecification. The companion RDS contains the
settings, fit-level parameter rows, and all summaries; no raw simulated data
are retained. The script measures Jacobian rank, condition number, loss,
convergence/admissibility, and alternative-basin flags rather than assuming the
constraint improves conditioning.

## Parameter-constraint results recorded on 2026-07-29

The command above was run with its displayed settings. Both exact signed-gamma
identities held to machine precision (maximum absolute error
`2.22e-16`). The two shape constraints reduced the independent parameter count
from 8 to 6 and increased nominal df from 4 to 6 in both kernels. All successful
fits had full column rank and finite SEs for every free and derived parameter.

| Kernel | Specification | Successful | Free | df | Rank | Median condition | Mean loss | Alternative-basin rate |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| contemporaneous | constrained | 50/50 | 6 | 6 | 6 | 72.5 | 0.0272 | 0.00 |
| contemporaneous | unconstrained | 50/50 | 8 | 4 | 8 | 16.2 | 0.0060 | 0.00 |
| dynamic | constrained | 48/50 | 6 | 6 | 6 | 191.6 | 0.0295 | 0.52 |
| dynamic | unconstrained | 50/50 | 8 | 4 | 8 | 109.8 | 0.0209 | 0.46 |

Thus the constraint increased df but did not improve conditioning in this
design. Two constrained dynamic replications produced no finite stationary
solution and were retained as non-admissible rather than omitted from the
convergence denominator.

For the contemporaneous constrained fits, bias/RMSE was `-0.0139/0.0520` and
`0.0080/0.0544` for the two directed paths. Shape RMSEs were 2.12 and 2.98;
derived skewness RMSEs were 0.200 and 0.198, and derived fourth-moment RMSEs
were 0.626 and 0.555. Empirical-SD/mean-SE ratios ranged from 1.30 to 1.80 and
coverage from 0.72 to 0.92. The Wilson interval for the lowest coverage was
`[0.583, 0.825]`.

Among the 48 admissible dynamic constrained fits, bias/RMSE for
`phi_X`, `X_lag_to_Y`, `Y_lag_to_X`, and `phi_Y` was respectively
`-0.0429/0.106`, `0.0293/0.142`, `-0.0409/0.181`, and `-0.0657/0.114`.
Shape RMSEs were 1.52 and 2.16. Derived third-cumulant RMSEs were 0.179 and
0.115; derived fourth-cumulant RMSEs were 0.573 and 0.271. Derived-parameter
empirical-SD/mean-SE ratios were 1.05--1.10, with coverage 0.81--0.88 and
Wilson intervals reported in the CSV. Transition-path sandwich SEs were
conservative in this modest design (ratios 0.46--0.59).

The deliberately wrong common-shape constraint raised contemporaneous loss to
0.744; the dynamic optimizer found no admissible stationary solution. The
strongly skewed non-gamma mixture raised contemporaneous loss to 0.349 and
produced large residual cumulant discrepancies; its dynamic fit was also
rejected by the stationarity/admissibility diagnostic. These are controlled
misspecification signals, not apparently exact recovery. With only 50
replications, coverage estimates remain noisy; the Wilson intervals should be
used instead of treating small deviations from 95% as definitive.
