# MCMSEM 0.27.0

## New stationary dynamic kernel

- Added `kernel = "dynamic"` to `MCMmodel()` for observed-state stationary
  VAR(1) moment models. `kernel = "contemporaneous"` remains the default and
  preserves the previous implementation numerically.
- Added `kernel = "static"` as a supported, silent alias that normalizes to
  `"contemporaneous"`; kernel values are otherwise matched exactly.
- Made the kernel choice explicit throughout the README, wiki, and package
  examples. Added reproducible longitudinal examples comparing CLPM/RI-CLPM
  estimates with dynamic MCMSEM fits to the largest complete cross-sectional
  wave.
- Added differentiable torch propagation of stationary cumulants through orders
  two, three, and four using linear solves with Kronecker powers of `B`.
- Fixed dynamic innovation variances to one and restricted innovation third and
  fourth cumulants to their diagonal entries.
- Added an optional full Gaussian residual/random-intercept covariance using a
  differentiable Cholesky parameterization. Gaussian residuals enter covariance
  and raw fourth-moment pairings, but not third or fourth cumulants.
- Added spectral-radius stationarity enforcement, nonnegative default bounds for
  autoregressions, signed cross-lagged paths, seeded multi-start fitting, and
  all-start diagnostics.
- Added dynamic fields to ordinary `mcmresultclass` objects, including `B`,
  `Psi_G`, innovation cumulants, spectral radius, stationarity, convergence,
  moment residuals, parameter/moment counts, and nominal degrees of freedom.
- Added `MCMmomentcount()`, `MCMdegreesoffreedom()`, `MCMdiagnostics()`, and
  `MCMdynamicloss()`.
- Added mathematical regression tests, gradient and permutation checks, a fast
  recovery simulation, an opt-in 20-replication recovery Monte Carlo script,
  and an opt-in 50-replication SE-calibration script.
- Added identity, diagonal, and full WLS moment weighting for dynamic fits.
  The prepared moment covariance now uses influence functions for estimated
  raw central moments through order four. Dynamic asymptotic covariance is
  available as the robust sandwich estimator or, with full WLS, the efficient
  inverse-information estimator. Derived `Psi_G` SEs use the delta method.
- Documentation credits both the original MCM-SEM framework to Tamimy, van
  Bergen, van der Zee, Dolan, and Nivard (2022), *Multi Co-Moment Structural
  Equation Models: Discovering Direction of Causality in the Presence of
  Confounding*, and the cumulant-identification/discrete-Lyapunov framework to
  Cecilie Olesen Recke, Sarah Lumpp, Nataliia Kushnerchuk, Janike Oldekop,
  Jiayi Li, Jane Ivy Coons, and Elina Robeva (2026), *Identifiability in
  Graphical Discrete Lyapunov Models*, arXiv preprint arXiv:2601.21818.

## Compatibility and fixes

- Existing model objects without kernel metadata are treated as
  contemporaneous models.
- Summary-data fits no longer require asymptotic-SE preparation when
  `compute_se = FALSE`.
- Corrected contemporaneous result extraction so `M4` is governed by
  `use_kurtosis`.
- Added the missing `methods` DESCRIPTION import and corrected the license-file
  spelling used by package checks.

## Known limitations

Dynamic WLS and asymptotic SEs currently require unweighted data and all moments
through order four. Older saved summaries must be rebuilt with MCMSEM 0.27.0 so
their moment covariance contains the required mean-estimation correction.
