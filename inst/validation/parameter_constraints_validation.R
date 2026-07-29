#!/usr/bin/env Rscript

# Opt-in validation for free, fixed, and derived MCMSEM parameters. This script
# uses synthetic signed-gamma and non-gamma data only. It deliberately does not
# read package examples, private data, or SIPP material.

parse_arguments <- function(args) {
  settings <- list(
    repetitions = 50L, n = 2000L, master_seed = 20260729L,
    rprop_iters = 80L, lbfgs_iters = 5L, starts = 2L,
    output = "validation-output/parameter-constraints"
  )
  for (argument in args) {
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !(pieces[1L] %in% names(settings))) {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
    name <- pieces[1L]
    value <- pieces[2L]
    settings[[name]] <- if (name == "output") value else as.integer(value)
  }
  settings
}

settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (settings$repetitions < 50L) {
  warning("Fewer than 50 replications were requested; use 50 or more for the recorded evaluation.")
}
dir.create(settings$output, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to validate the MCMSEM working tree.")
}
devtools::load_all(".", quiet = TRUE)

signed_gamma <- function(n, shape, sign) {
  sign * (stats::rgamma(n, shape = shape) - shape) / sqrt(shape)
}

non_gamma_innovation <- function(n, sign) {
  # A strongly skewed Gaussian mixture whose third/fourth cumulants do not
  # satisfy kappa_4 = 1.5 * kappa_3^2 (the gap is about -2.6 here).
  component <- stats::rbinom(n, 1, 0.10)
  x <- stats::rnorm(n, mean = 3 * component, sd = 0.4 + 0.4 * component)
  sign * as.numeric(scale(x))
}

simulate_contemporaneous <- function(n, A, shapes, signs,
                                     family = c("gamma", "mixture")) {
  family <- match.arg(family)
  innovations <- vapply(seq_along(shapes), function(i) {
    if (family == "gamma") signed_gamma(n, shapes[i], signs[i]) else
      non_gamma_innovation(n, signs[i])
  }, numeric(n))
  observed <- innovations %*% t(solve(diag(nrow(A)) - A))
  as.data.frame(observed)
}

simulate_dynamic <- function(n, B, shapes, signs,
                             family = c("gamma", "mixture"), burnin = 180L) {
  family <- match.arg(family)
  state <- matrix(0, n, nrow(B))
  for (time in seq_len(burnin)) {
    innovations <- vapply(seq_along(shapes), function(i) {
      if (family == "gamma") signed_gamma(n, shapes[i], signs[i]) else
        non_gamma_innovation(n, signs[i])
    }, numeric(n))
    state <- state %*% t(B) + innovations
  }
  as.data.frame(state)
}

set_starts <- function(model, values) {
  for (name in intersect(names(values), model$param_names)) {
    model <- MCMedit(model, "start", name, values[[name]])
  }
  model
}

add_gamma_constraints <- function(model, kernel, shapes = c(X = 4, Y = 6),
                                  signs = c(X = -1, Y = 1), common = FALSE) {
  if (common) {
    model <- MCMparameter(
      model, "shape_common", "free", start = mean(shapes),
      transform = "positive"
    )
  }
  for (variable in names(shapes)) {
    shape_name <- if (common) "shape_common" else paste0("shape_", variable)
    if (!common) {
      model <- MCMparameter(
        model, shape_name, "free", start = shapes[[variable]],
        transform = "positive"
      )
    }
    sign_name <- paste0("sign_", variable)
    model <- MCMparameter(model, sign_name, "fixed", value = signs[[variable]])
  }
  if (kernel == "dynamic") {
    model <- MCMparameter(
      model, "tau_X", "derived",
      expression = if (common) {
        ~ sign_X * 2 / sqrt(shape_common)
      } else ~ sign_X * 2 / sqrt(shape_X)
    )
    model <- MCMparameter(
      model, "tau_Y", "derived",
      expression = if (common) {
        ~ sign_Y * 2 / sqrt(shape_common)
      } else ~ sign_Y * 2 / sqrt(shape_Y)
    )
    model <- MCMparameter(
      model, "kappa_X", "derived",
      expression = if (common) ~ 6 / shape_common else ~ 6 / shape_X
    )
    model <- MCMparameter(
      model, "kappa_Y", "derived",
      expression = if (common) ~ 6 / shape_common else ~ 6 / shape_Y
    )
  } else {
    model <- MCMparameter(
      model, "sk1", "derived",
      expression = if (common) {
        ~ sign_X * 2 / sqrt(shape_common)
      } else ~ sign_X * 2 / sqrt(shape_X)
    )
    model <- MCMparameter(
      model, "sk2", "derived",
      expression = if (common) {
        ~ sign_Y * 2 / sqrt(shape_common)
      } else ~ sign_Y * 2 / sqrt(shape_Y)
    )
    model <- MCMparameter(
      model, "k1", "derived",
      expression = if (common) ~ 3 + 6 / shape_common else ~ 3 + 6 / shape_X
    )
    model <- MCMparameter(
      model, "k2", "derived",
      expression = if (common) ~ 3 + 6 / shape_common else ~ 3 + 6 / shape_Y
    )
  }
  model
}

model_set <- function(summary_data, kernel, shapes = c(X = 4, Y = 6),
                      signs = c(X = -1, Y = 1)) {
  if (kernel == "dynamic") {
    base <- MCMmodel(summary_data, kernel = "dynamic", gaussian_residual = FALSE)
    truth <- c(
      phi_X = 0.42, X_lag_to_Y = -0.11,
      Y_lag_to_X = 0.16, phi_Y = 0.32,
      tau_X = signs[["X"]] * 2 / sqrt(shapes[["X"]]),
      tau_Y = signs[["Y"]] * 2 / sqrt(shapes[["Y"]]),
      kappa_X = 6 / shapes[["X"]], kappa_Y = 6 / shapes[["Y"]]
    )
  } else {
    base <- MCMmodel(summary_data, n_latent = 0, kernel = "contemporaneous")
    truth <- c(
      b1_2 = -0.10, b2_1 = 0.15, s1 = 1, s2 = 1,
      sk1 = signs[["X"]] * 2 / sqrt(shapes[["X"]]),
      sk2 = signs[["Y"]] * 2 / sqrt(shapes[["Y"]]),
      k1 = 3 + 6 / shapes[["X"]],
      k2 = 3 + 6 / shapes[["Y"]]
    )
  }
  unconstrained <- set_starts(base, truth)
  constrained <- add_gamma_constraints(base, kernel, shapes, signs, common = FALSE)
  constrained <- set_starts(constrained, truth)
  misspecified <- add_gamma_constraints(base, kernel, shapes, signs, common = TRUE)
  misspecified <- set_starts(misspecified, truth)
  list(
    unconstrained = unconstrained, constrained = constrained,
    misspecified = misspecified, truth = truth
  )
}

simulate_from_model <- function(model, n, shapes, signs,
                                family = c("gamma", "mixture")) {
  family <- match.arg(family)
  if (.model_kernel(model) == "dynamic") {
    simulate_dynamic(n, model$num_matrices$B, shapes, signs, family)
  } else {
    simulate_contemporaneous(n, model$num_matrices$A, shapes, signs, family)
  }
}

safe_fit <- function(model, data, kernel, compute_se = TRUE) {
  # MCMfit seeds multi-start generation. Preserve the outer Monte Carlo stream
  # so fitting one replication cannot determine the data generated in the next.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) simulation_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", simulation_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  tryCatch(
    do.call(MCMfit, c(list(
      mcmmodel = model, data = data, compute_se = compute_se,
      optimizers = c("rprop", "lbfgs"),
      optim_iters = c(settings$rprop_iters, settings$lbfgs_iters),
      learning_rate = c(0.01, 0.3), n_starts = settings$starts,
      seed = settings$master_seed, moment_weighting = "identity"
    ), if (kernel == "dynamic") list(se_correction = "robust") else list())),
    error = function(e) structure(
      list(error = conditionMessage(e)), class = "constraint_validation_error"
    )
  )
}

fit_diagnostics <- function(fit) {
  if (inherits(fit, "constraint_validation_error")) {
    return(list(rank = NA_integer_, condition = NA_real_, stationary = FALSE,
                alternative_basin = NA))
  }
  diagnostics <- tryCatch(MCMdiagnostics(fit, jacobian = TRUE),
                          error = function(e) NULL)
  alternative_basin <- if (.result_kernel(fit) == "dynamic" &&
                            nrow(fit$start_diagnostics) > 1L) {
    admissible <- fit$start_diagnostics$final_loss[
      fit$start_diagnostics$convergence == 0L
    ]
    length(admissible) > 1L && diff(range(admissible)) > 1e-6
  } else FALSE
  list(
    rank = if (is.null(diagnostics)) NA_integer_ else diagnostics$jacobian_rank,
    condition = if (is.null(diagnostics)) NA_real_ else diagnostics$jacobian_condition,
    stationary = if (.result_kernel(fit) == "dynamic") isTRUE(fit$stationary) else TRUE,
    alternative_basin = alternative_basin
  )
}

truth_for_model <- function(model, shapes, signs) {
  values <- stats::setNames(MCMparameters(model)$value, MCMparameters(model)$name)
  truth <- values
  if (.model_kernel(model) == "dynamic") {
    truth[c("phi_X", "X_lag_to_Y", "Y_lag_to_X", "phi_Y")] <-
      c(0.42, -0.11, 0.16, 0.32)
    truth[c("tau_X", "tau_Y")] <- signs * 2 / sqrt(shapes)
    truth[c("kappa_X", "kappa_Y")] <- 6 / shapes
  } else {
    truth[c("b1_2", "b2_1", "s1", "s2")] <- c(-0.10, 0.15, 1, 1)
    truth[c("sk1", "sk2")] <- signs * 2 / sqrt(shapes)
    truth[c("k1", "k2")] <- 3 + 6 / shapes
  }
  for (name in intersect(c("shape_X", "shape_Y"), names(truth))) {
    truth[[name]] <- shapes[[sub("shape_", "", name)]]
  }
  truth
}

fit_rows <- function(fit, kernel, specification, scenario, replication,
                     shapes, signs) {
  if (inherits(fit, "constraint_validation_error")) {
    return(data.frame(
      kernel = kernel, specification = specification, scenario = scenario,
      replication = replication, parameter = NA_character_, type = NA_character_,
      truth = NA_real_, estimate = NA_real_, se = NA_real_, covered = NA,
      loss = NA_real_, n_free = NA_integer_, df = NA_integer_, rank = NA_integer_,
      condition = NA_real_, converged = FALSE, admissible = FALSE,
      alternative_basin = NA, error = fit$error, stringsAsFactors = FALSE
    ))
  }
  table <- fit$parameter_table
  truth <- truth_for_model(fit$model, shapes, signs)
  diagnostics <- fit_diagnostics(fit)
  dof <- MCMdegreesoffreedom(fit)
  se <- table$se
  covered <- is.finite(se) &
    table$estimate - 1.96 * se <= truth[table$parameter] &
    table$estimate + 1.96 * se >= truth[table$parameter]
  data.frame(
    kernel = kernel, specification = specification, scenario = scenario,
    replication = replication, parameter = table$parameter, type = table$type,
    truth = unname(truth[table$parameter]), estimate = table$estimate,
    se = se, covered = covered, loss = fit$loss,
    n_free = dof$n_parameters, df = dof$df, rank = diagnostics$rank,
    condition = diagnostics$condition, converged = is.finite(fit$loss),
    admissible = diagnostics$stationary, alternative_basin = diagnostics$alternative_basin,
    error = NA_character_, stringsAsFactors = FALSE
  )
}

summarize_monte_carlo <- function(rows) {
  rows <- rows[rows$type %in% c("free", "derived") & rows$converged, ]
  groups <- split(rows, interaction(rows$kernel, rows$specification,
                                    rows$scenario, rows$parameter, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    estimate <- x$estimate
    se <- x$se
    coverage <- mean(x$covered, na.rm = TRUE)
    n_coverage <- sum(!is.na(x$covered))
    z <- 1.96
    denominator <- 1 + z^2 / n_coverage
    center <- (coverage + z^2 / (2 * n_coverage)) / denominator
    half <- z * sqrt(coverage * (1 - coverage) / n_coverage +
                       z^2 / (4 * n_coverage^2)) / denominator
    data.frame(
      kernel = x$kernel[1], specification = x$specification[1],
      scenario = x$scenario[1], parameter = x$parameter[1], type = x$type[1],
      true = x$truth[1], mean_estimate = mean(estimate),
      empirical_sd = stats::sd(estimate), mean_se = mean(se, na.rm = TRUE),
      sd_to_mean_se = stats::sd(estimate) / mean(se, na.rm = TRUE),
      bias = mean(estimate - x$truth),
      rmse = sqrt(mean((estimate - x$truth)^2)), coverage = coverage,
      coverage_wilson_lower = max(0, center - half),
      coverage_wilson_upper = min(1, center + half),
      replications = length(unique(x$replication)), stringsAsFactors = FALSE
    )
  }))
}

summarize_specification <- function(rows) {
  unique_fit <- rows[!duplicated(rows[c("kernel", "specification", "scenario",
                                        "replication")]), ]
  groups <- split(unique_fit, interaction(unique_fit$kernel,
                                           unique_fit$specification,
                                           unique_fit$scenario, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) data.frame(
    kernel = x$kernel[1], specification = x$specification[1],
    scenario = x$scenario[1], n_free = mean(x$n_free, na.rm = TRUE),
    df = mean(x$df, na.rm = TRUE),
    jacobian_rank = mean(x$rank, na.rm = TRUE),
    median_condition = stats::median(x$condition, na.rm = TRUE),
    mean_loss = mean(x$loss, na.rm = TRUE),
    alternative_basin_rate = mean(x$alternative_basin, na.rm = TRUE),
    convergence_rate = mean(x$converged),
    admissibility_rate = mean(x$admissible), stringsAsFactors = FALSE
  )))
}

shapes <- c(X = 4, Y = 6)
signs <- c(X = -1, Y = 1)
exact <- do.call(rbind, lapply(names(shapes), function(variable) {
  tau <- signs[[variable]] * 2 / sqrt(shapes[[variable]])
  kappa <- 6 / shapes[[variable]]
  data.frame(
    variable = variable, shape = shapes[[variable]], sign = signs[[variable]],
    tau = tau, kappa = kappa, gamma_identity = 1.5 * tau^2,
    absolute_error = abs(kappa - 1.5 * tau^2), stringsAsFactors = FALSE
  )
}))
stopifnot(max(exact$absolute_error) < 1e-10)
utils::write.csv(exact, file.path(settings$output, "exact_moment_checks.csv"),
                 row.names = FALSE)

set.seed(settings$master_seed)
all_rows <- list()
row_index <- 1L
for (kernel in c("contemporaneous", "dynamic")) {
  message("Validating ", kernel, " kernel")
  template_data <- as.data.frame(matrix(stats::rnorm(4000), ncol = 2))
  names(template_data) <- names(shapes)
  template_summary <- MCMdatasummary(
    template_data, scale_data = FALSE, prep_asymptotic_se = TRUE
  )
  templates <- model_set(template_summary, kernel, shapes, signs)
  generating_model <- templates$constrained

  for (replication in seq_len(settings$repetitions)) {
    if (replication %% 5L == 0L) {
      message("  replication ", replication, "/", settings$repetitions)
    }
    raw <- simulate_from_model(
      generating_model, settings$n, shapes, signs, family = "gamma"
    )
    names(raw) <- names(shapes)
    summary_data <- MCMdatasummary(
      raw, scale_data = FALSE, prep_asymptotic_se = TRUE
    )
    models <- model_set(summary_data, kernel, shapes, signs)
    for (specification in c("constrained", "unconstrained")) {
      fit <- safe_fit(models[[specification]], summary_data, kernel, TRUE)
      all_rows[[row_index]] <- fit_rows(
        fit, kernel, specification, "gamma", replication, shapes, signs
      )
      row_index <- row_index + 1L
    }
  }

  # Deliberately wrong common-shape constraint on correctly generated data.
  raw <- simulate_from_model(
    generating_model, settings$n, shapes, signs, family = "gamma"
  )
  names(raw) <- names(shapes)
  summary_data <- MCMdatasummary(raw, scale_data = FALSE,
                                 prep_asymptotic_se = TRUE)
  misspecified_models <- model_set(summary_data, kernel, shapes, signs)
  fit <- safe_fit(misspecified_models$misspecified, summary_data, kernel, TRUE)
  all_rows[[row_index]] <- fit_rows(
    fit, kernel, "common_shape", "wrong_constraint", 1L, shapes, signs
  )
  row_index <- row_index + 1L

  # Stable framework behavior under non-gamma innovation misspecification.
  raw <- simulate_from_model(
    generating_model, settings$n, shapes, signs, family = "mixture"
  )
  names(raw) <- names(shapes)
  summary_data <- MCMdatasummary(raw, scale_data = FALSE,
                                 prep_asymptotic_se = TRUE)
  mixture_models <- model_set(summary_data, kernel, shapes, signs)
  fit <- safe_fit(mixture_models$constrained, summary_data, kernel, TRUE)
  all_rows[[row_index]] <- fit_rows(
    fit, kernel, "constrained", "non_gamma", 1L, shapes, signs
  )
  row_index <- row_index + 1L
}

rows <- do.call(rbind, all_rows)
monte_carlo <- summarize_monte_carlo(rows[rows$scenario == "gamma", ])
constraint_benefit <- summarize_specification(rows[rows$scenario == "gamma", ])
recovery <- aggregate(
  cbind(bias = rows$estimate - rows$truth,
        squared_error = (rows$estimate - rows$truth)^2) ~
    kernel + specification + scenario + parameter + type,
  data = rows, FUN = function(x) mean(x, na.rm = TRUE)
)
recovery$rmse <- sqrt(recovery$squared_error)
misspecification <- summarize_specification(
  rows[rows$scenario %in% c("wrong_constraint", "non_gamma"), ]
)

utils::write.csv(recovery, file.path(settings$output, "recovery_summary.csv"),
                 row.names = FALSE)
utils::write.csv(monte_carlo,
                 file.path(settings$output, "monte_carlo_se_summary.csv"),
                 row.names = FALSE)
utils::write.csv(constraint_benefit,
                 file.path(settings$output, "constraint_benefit_summary.csv"),
                 row.names = FALSE)
utils::write.csv(misspecification,
                 file.path(settings$output, "misspecification_summary.csv"),
                 row.names = FALSE)
saveRDS(
  list(settings = settings, exact = exact, fit_rows = rows,
       recovery = recovery, monte_carlo = monte_carlo,
       constraint_benefit = constraint_benefit,
       misspecification = misspecification),
  file.path(settings$output, "parameter_constraint_validation.rds")
)

message("Validation outputs written to ", normalizePath(settings$output))
