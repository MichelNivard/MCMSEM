# Opt-in methodological validation

These scripts are opt-in methodological checks and are not run during package
tests. Generated `validation-output/` files are deliberately excluded from
source builds and Git.

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
