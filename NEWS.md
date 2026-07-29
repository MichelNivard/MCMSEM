# MCMSEM 0.27.0

## General parameter constraints

- Added `MCMparameter()` and `MCMparameters()` with one canonical graph for
  independent free parameters, fixed constants, and safely derived parameters.
  Auxiliary free parameters need not occupy matrix cells, and derived
  parameters can be chained or reused across cells in either kernel.
- Added a restricted differentiable expression language supporting symbols,
  finite numeric constants, arithmetic, powers, `sqrt()`, `exp()`, `log()`,
  `softplus()`, and `logistic()`. Expressions are compiled separately to base R
  and Torch without unrestricted parsing or evaluation; unknown symbols,
  unsupported calls, dependency cycles, invalid domains, and non-finite starts
  are rejected.
- Added exact positive (softplus) and finite-interval (logistic)
  transformations. Natural reported values remain separate from internal
  optimizer coordinates. Existing parameters retain identity coordinates and
  their historical bound-penalty behavior.
- Only independent free parameters now enter optimization, parameter counts,
  degrees of freedom, Jacobians, gradient histories, and information criteria.
  Fixed and derived parameters remain in result tables and summaries without
  consuming degrees of freedom.
- Added delta-method covariance and SE propagation from optimizer coordinates
  to natural free and derived parameters. Fixed-parameter SEs are zero when an
  estimator covariance is available; all parameter covariances are retained on
  fitted result objects.
- Added `MCMimpliedmoments()` as a common base-R evaluator for diagnostics and
  validation in both kernels. Torch/base parity and gradient/finite-difference
  tests cover the same graph compiler.
- Extended copying, RDS serialization, result refitting, diagnostics, model
  comparison, and summaries to retain and report parameter types,
  transformations, dependencies, expressions, and matrix locations.
- Added signed-gamma examples with positive shape, fixed sign, and derived
  third/fourth cumulants. In the contemporaneous kernel, `K` stores the raw
  standardized fourth moment (`3 + 6 / shape`); dynamic `Kappa` stores the
  fourth cumulant (`6 / shape`).

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
- Fixed one-step bootstrap fitting with data frames by passing numeric matrices
  to the higher-order-moment routines; bootstrap covariance now propagates to
  free, fixed, and derived result parameters.
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
