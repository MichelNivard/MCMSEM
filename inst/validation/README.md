# Dynamic-kernel validation

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
