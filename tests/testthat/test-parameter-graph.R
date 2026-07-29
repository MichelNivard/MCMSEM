test_that("free, fixed, derived, chained, and auxiliary parameters are canonical", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(
    model, "shape_A", "free", start = 4, transform = "positive"
  )
  model <- MCMparameter(model, "sign_A", "fixed", value = -1)
  model <- MCMparameter(
    model, "tau_A", "derived", expression = ~ sign_A * 2 / sqrt(shape_A)
  )
  model <- MCMparameter(
    model, "kappa_A", "derived", expression = ~ 1.5 * tau_A^2
  )

  parameters <- MCMparameters(model)
  expect_identical(parameters$type[match(
    c("shape_A", "sign_A", "tau_A", "kappa_A"), parameters$name
  )], c("free", "fixed", "derived", "derived"))
  expect_equal(parameters$value[match("tau_A", parameters$name)], -1)
  expect_equal(parameters$value[match("kappa_A", parameters$name)], 1.5)
  expect_true(parameters$auxiliary[match("shape_A", parameters$name)])
  expect_false("shape_A" %in% unlist(model$named_matrices, use.names = FALSE))
  expect_identical(model$free_param_names, model$param_names)
  expect_equal(model$num_matrices$Tau[1, 1], -1)
  expect_equal(model$num_matrices$Kappa[1, 1], 1.5)
})

test_that("the safe expression language implements every supported operation", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "x", "free", start = 4)
  model <- MCMparameter(model, "constant", "fixed", value = 2)
  definitions <- list(
    plus = ~ x + constant,
    unary_plus = ~ +x,
    minus = ~ x - constant,
    unary_minus = ~ -x,
    multiply = ~ x * constant,
    divide = ~ x / constant,
    power = ~ x^constant,
    square_root = ~ sqrt(x),
    exponential = ~ exp(constant),
    logarithm = ~ log(x),
    soft_plus = ~ softplus(-x),
    logistic_value = ~ logistic(-x),
    numeric_constant = ~ 6 / 3
  )
  for (name in names(definitions)) {
    model <- MCMparameter(
      model, name, "derived", expression = definitions[[name]]
    )
  }
  values <- stats::setNames(MCMparameters(model)$value, MCMparameters(model)$name)
  expect_equal(values[c("plus", "unary_plus", "minus", "unary_minus")],
               c(plus = 6, unary_plus = 4, minus = 2, unary_minus = -4))
  expect_equal(values[c("multiply", "divide", "power", "square_root")],
               c(multiply = 8, divide = 2, power = 16, square_root = 2))
  expect_equal(values[["exponential"]], exp(2))
  expect_equal(values[["logarithm"]], log(4))
  expect_equal(values[["soft_plus"]], log1p(exp(-4)))
  expect_equal(values[["logistic_value"]], plogis(-4))
  expect_equal(values[["numeric_constant"]], 2)
})

test_that("invalid and conflicting graph definitions are rejected", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "x", "free", start = 1)
  expect_error(MCMparameter(model, "x", "free", start = 2), "already defined")
  expect_error(
    MCMparameter(model, "bad", "derived", expression = ~ system(x)),
    "Unsupported"
  )
  expect_error(
    MCMparameter(model, "bad", "derived", expression = ~ unknown + x),
    "Unknown symbol"
  )
  negative <- MCMparameter(model, "negative", "fixed", value = -1)
  expect_error(
    MCMparameter(negative, "bad", "derived", expression = ~ log(negative)),
    "non-finite"
  )
  cyclic <- MCMparameter(model, "y", "free", start = 2)
  cyclic <- MCMparameter(cyclic, "x", "derived", expression = ~ y + 1)
  expect_error(
    MCMparameter(cyclic, "y", "derived", expression = ~ x + 1),
    "cycle"
  )
  expect_error(
    MCMparameter(model, "bounded", "free", start = 1,
                 transform = "bounded", lower = 0, upper = 1),
    "strictly between"
  )
  expect_error(
    MCMparameter(model, "positive", "free", start = 0,
                 transform = "positive"),
    "above its lower bound"
  )
})

test_that("free parameters use invertible natural-to-optimizer transformations", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(
    model, "positive", "free", start = 0.25, transform = "positive"
  )
  model <- MCMparameter(
    model, "probability", "free", start = 0.2, transform = "bounded",
    lower = 0, upper = 1
  )
  eta <- MCMSEM:::.parameter_optimizer_coordinates(model)
  reported <- MCMSEM:::.parameter_values_base(
    model, eta, optimizer_scale = TRUE
  )
  expect_equal(reported[c("positive", "probability")],
               c(positive = 0.25, probability = 0.2), tolerance = 1e-12)
  updated <- model$param_values
  updated[match(c("positive", "probability"), model$param_names)] <- c(0.5, 0.3)
  model$param_values <- updated
  model$inverse_parse()
  expect_equal(
    model$parameter_table$optimizer_value[model$parameter_table$type == "free"],
    unname(MCMSEM:::.parameter_optimizer_coordinates(model)),
    tolerance = 1e-12
  )
})

test_that("parameters convert among types and invalid legacy edits are informative", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "tau_A", "fixed", value = -1)
  expect_false("tau_A" %in% model$param_names)
  expect_error(MCMedit(model, "start", "tau_A", 0), "only be edited for free")
  expect_error(MCMedit(model, "lbound", "tau_A", -2), "only be edited for free")
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(
    model, "tau_A", "derived", expression = ~ -2 / sqrt(shape)
  )
  expect_equal(model$num_matrices$Tau[1, 1], -1)
  model <- MCMparameter(model, "tau_A", "free", start = -0.8)
  expect_true("tau_A" %in% model$param_names)
  expect_equal(model$num_matrices$Tau[1, 1], -0.8)
})

test_that("MCMedit supports arbitrary equality labels with unbounded starts", {
  model <- MCMmodel(make_moment_summary())
  original <- model$named_matrices$A[1, 2]
  model <- MCMedit(model, "A", original, "custom_path")
  expect_true("custom_path" %in% model$param_names)
  expect_true(MCMSEM:::.parameter_graph_active(model))
  bounds <- MCMSEM:::.parameter_optimizer_bounds(model)
  expect_true(all(is.finite(c(bounds$L, bounds$U))))
  expect_equal(model$num_matrices$A[1, 2], 0)
})

test_that("fully fixed models fit without optimizer coordinates", {
  for (kernel in c("contemporaneous", "dynamic")) {
    model <- MCMmodel(make_moment_summary(), kernel = kernel)
    for (name in model$param_names) {
      model <- MCMparameter(
        model, name, "fixed",
        value = model$param_values[match(name, model$param_names)]
      )
    }
    expect_length(model$param_names, 0L)
    fit <- MCMfit(model, make_moment_summary(), compute_se = FALSE)
    expect_length(fit$model$param_names, 0L)
    expect_true(is.finite(fit$loss))
    expect_true(all(fit$parameter_table$type == "fixed"))
    diagnostics <- MCMdiagnostics(model, jacobian = TRUE)
    expect_equal(diagnostics$jacobian_rank, 0L)
    expect_identical(dim(diagnostics$jacobian), c(12L, 0L))
    expect_true(is.na(diagnostics$jacobian_condition))
  }
})

test_that("derived equality labels remain one reported quantity", {
  model <- MCMmodel(make_moment_summary())
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sk1", "derived",
                        expression = ~ 2 / sqrt(shape))
  model <- MCMedit(model, "Sk", "sk2", "sk1")
  row <- match("sk1", model$parameter_table$name)
  expect_equal(nrow(model$parameter_table$matrix_locations[[row]]), 2L)
  expect_equal(sum(model$parameter_table$name == "sk1"), 1L)
  expect_equal(model$num_matrices$Sk[
    is.na(suppressWarnings(as.numeric(model$named_matrices$Sk)))
  ], c(1, 1))
})

test_that("copying and serialization retain an independent parameter graph", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ -2 / sqrt(shape))
  copied <- model$copy()
  copied <- MCMedit(copied, "start", "shape", 9)
  expect_equal(MCMparameters(model)$value[MCMparameters(model)$name == "shape"], 4)
  expect_equal(MCMparameters(copied)$value[MCMparameters(copied)$name == "shape"], 9)

  path <- tempfile(fileext = ".rds")
  saveRDS(model, path)
  restored <- readRDS(path)
  expect_equal(MCMparameters(restored), MCMparameters(model))
  expect_equal(MCMimpliedmoments(restored)$M3, MCMimpliedmoments(model)$M3)
})

test_that("signed-gamma constraints reduce only the independent free count", {
  unconstrained <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  constrained <- unconstrained
  for (variable in c("A", "B")) {
    shape <- paste0("shape_", variable)
    sign <- paste0("sign_", variable)
    constrained <- MCMparameter(
      constrained, shape, "free", start = if (variable == "A") 4 else 6,
      transform = "positive"
    )
    constrained <- MCMparameter(
      constrained, sign, "fixed", value = if (variable == "A") -1 else 1
    )
    constrained <- MCMparameter(
      constrained, paste0("tau_", variable), "derived",
      expression = if (variable == "A") {
        ~ sign_A * 2 / sqrt(shape_A)
      } else {
        ~ sign_B * 2 / sqrt(shape_B)
      }
    )
    constrained <- MCMparameter(
      constrained, paste0("kappa_", variable), "derived",
      expression = if (variable == "A") ~ 6 / shape_A else ~ 6 / shape_B
    )
  }
  old_df <- MCMdegreesoffreedom(unconstrained)
  new_df <- MCMdegreesoffreedom(constrained)
  expect_equal(old_df$n_parameters - new_df$n_parameters, 2)
  expect_equal(new_df$df - old_df$df, 2)
  values <- stats::setNames(MCMparameters(constrained)$value,
                           MCMparameters(constrained)$name)
  expect_equal(values[["kappa_A"]], 1.5 * values[["tau_A"]]^2,
               tolerance = 1e-12)
  expect_equal(values[["kappa_B"]], 1.5 * values[["tau_B"]]^2,
               tolerance = 1e-12)
})

test_that("base R and Torch graph evaluation agree in both kernels", {
  dynamic <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  dynamic <- MCMparameter(dynamic, "shape_A", "free", start = 4,
                          transform = "positive")
  dynamic <- MCMparameter(dynamic, "sign_A", "fixed", value = -1)
  dynamic <- MCMparameter(dynamic, "tau_A", "derived",
                          expression = ~ sign_A * 2 / sqrt(shape_A))
  dynamic <- MCMparameter(dynamic, "kappa_A", "derived",
                          expression = ~ 6 / shape_A)
  components <- MCMSEM:::.get_dynamic_torch_components(
    dynamic, torch_device("cpu"), torch_float64()
  )
  torch_implied <- MCMSEM:::.get_dynamic_predicted_matrices(components)
  base_implied <- MCMimpliedmoments(dynamic)
  expect_equal(as.matrix(torch_implied$M2), base_implied$M2, tolerance = 1e-12)
  expect_equal(as.matrix(torch_implied$M3), base_implied$M3, tolerance = 1e-12)
  expect_equal(as.matrix(torch_implied$M4), base_implied$M4, tolerance = 1e-11)

  contemporaneous <- MCMmodel(make_moment_summary())
  contemporaneous <- MCMparameter(
    contemporaneous, "shape", "free", start = 4, transform = "positive"
  )
  contemporaneous <- MCMparameter(contemporaneous, "sign", "fixed", value = -1)
  contemporaneous <- MCMparameter(
    contemporaneous, "sk1", "derived", expression = ~ sign * 2 / sqrt(shape)
  )
  contemporaneous <- MCMparameter(
    contemporaneous, "k1", "derived", expression = ~ 3 + 6 / shape
  )
  observed <- list(
    M2 = torch_tensor(matrix(0, 2, 2)),
    M3 = torch_tensor(matrix(0, 2, 4)),
    M4 = torch_tensor(matrix(0, 2, 8))
  )
  matrices <- MCMSEM:::.get_torch_matrices(
    contemporaneous, torch_device("cpu"), observed$M2, observed$M3,
    observed$M4, torch_float64()
  )
  torch_contemporaneous <- MCMSEM:::.get_predicted_matrices(
    matrices$.par_list, matrices$torch_masks, matrices$torch_maps,
    matrices$base_matrices, TRUE, TRUE, TRUE, FALSE,
    jit_compile(MCMSEM:::.jit_funcs$slownecker)
  )
  base_contemporaneous <- MCMimpliedmoments(contemporaneous)
  expect_equal(as.matrix(torch_contemporaneous$M2),
               base_contemporaneous$M2, tolerance = 1e-12)
  expect_equal(as.matrix(torch_contemporaneous$M3),
               base_contemporaneous$M3, tolerance = 1e-12)
  expect_equal(as.matrix(torch_contemporaneous$M4),
               base_contemporaneous$M4, tolerance = 1e-11)
})

test_that("Torch graph gradients agree with finite differences", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "shape_A", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign_A", "fixed", value = -1)
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ sign_A * 2 / sqrt(shape_A))
  model <- MCMparameter(model, "kappa_A", "derived",
                        expression = ~ 6 / shape_A)
  eta <- MCMSEM:::.parameter_optimizer_coordinates(model)
  tensor <- torch_tensor(unname(eta), requires_grad = TRUE, dtype = torch_float64())
  values <- MCMSEM:::.parameter_values_torch(model, tensor)
  objective <- values$tau_A + values$kappa_A
  objective$backward()
  torch_gradient <- as.numeric(tensor$grad)
  finite_difference <- numDeriv::grad(function(x) {
    values <- MCMSEM:::.parameter_values_base(model, x, optimizer_scale = TRUE)
    values[["tau_A"]] + values[["kappa_A"]]
  }, unname(eta))
  expect_equal(torch_gradient, finite_difference, tolerance = 1e-6)
})

test_that("delta-method covariance matches analytic signed-gamma derivatives", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "shape_A", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign_A", "fixed", value = -1)
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ sign_A * 2 / sqrt(shape_A))
  model <- MCMparameter(model, "kappa_A", "derived",
                        expression = ~ 6 / shape_A)
  eta <- MCMSEM:::.parameter_optimizer_coordinates(model)
  V_eta <- diag(length(eta)) * 0
  dimnames(V_eta) <- list(model$param_names, model$param_names)
  shape_index <- match("shape_A", model$param_names)
  V_eta[shape_index, shape_index] <- 0.04
  delta <- MCMSEM:::.parameter_covariance_from_optimizer(
    model, V_eta, eta, method = "Richardson"
  )
  d_shape <- plogis(eta[["shape_A"]])
  d_tau <- d_shape / 4^(3 / 2)
  d_kappa <- -6 * d_shape / 4^2
  expect_equal(delta$vcov["shape_A", "shape_A"], d_shape^2 * 0.04,
               tolerance = 1e-8)
  expect_equal(delta$vcov["tau_A", "tau_A"], d_tau^2 * 0.04,
               tolerance = 1e-8)
  expect_equal(delta$vcov["kappa_A", "kappa_A"], d_kappa^2 * 0.04,
               tolerance = 1e-8)
  expect_equal(delta$se[["sign_A"]], 0)
})

test_that("contemporaneous asymptotic covariance propagates to reported parameters", {
  model <- MCMmodel(make_moment_summary())
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign", "fixed", value = -1)
  model <- MCMparameter(model, "sk1", "derived",
                        expression = ~ sign * 2 / sqrt(shape))
  model <- MCMparameter(model, "k1", "derived",
                        expression = ~ 3 + 6 / shape)
  data <- prepare_population_vcov(make_moment_summary(), diag(12) / 10000)
  observed <- list(
    M2 = torch_tensor(data$M2, dtype = torch_float32()),
    M3 = torch_tensor(data$M3, dtype = torch_float32()),
    M4 = torch_tensor(data$M4, dtype = torch_float32())
  )
  matrices <- MCMSEM:::.get_torch_matrices(
    model, torch_device("cpu"), observed$M2, observed$M3, observed$M4,
    torch_float32()
  )
  coordinate_covariance <- MCMSEM:::.std.err(
    data, matrices$.par_list, TRUE, TRUE, matrices$torch_masks,
    matrices$torch_maps, matrices$base_matrices, matrices$m2v_masks,
    torch_device("cpu"), FALSE, TRUE, "simple", FALSE,
    return_vcov = TRUE
  )$vcov
  reported <- MCMSEM:::.parameter_covariance_from_optimizer(
    model, coordinate_covariance,
    MCMSEM:::.parameter_optimizer_coordinates(model)
  )
  expect_true(all(is.finite(reported$se)))
  expect_equal(reported$se[["sign"]], 0)
  expect_gt(reported$se[["sk1"]], 0)
  expect_gt(reported$se[["k1"]], 0)
})

test_that("contemporaneous bootstrap covariance includes derived parameters", {
  set.seed(20260729)
  raw <- data.frame(A = rnorm(1200), B = rnorm(1200))
  summary_data <- suppressWarnings(MCMdatasummary(
    raw, prep_asymptotic_se = FALSE
  ))
  model <- MCMmodel(summary_data)
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign", "fixed", value = -1)
  model <- MCMparameter(model, "sk1", "derived",
                        expression = ~ sign * 2 / sqrt(shape))
  model <- MCMparameter(model, "k1", "derived",
                        expression = ~ 3 + 6 / shape)
  for (se_type in c("one-step", "two-step")) {
    invisible(capture.output(fit <- suppressWarnings(MCMfit(
      model, raw, compute_se = TRUE, se_type = se_type,
      bootstrap_iter = 3, bootstrap_chunks = 20,
      optimizers = "rprop", optim_iters = 1, learning_rate = 0.01
    ))))
    expect_identical(dim(fit$parameter_vcov),
                     rep(nrow(fit$parameter_table), 2L))
    expect_equal(fit$parameter_table$se[
      fit$parameter_table$parameter == "sign"
    ], 0)
    expect_true(all(is.finite(fit$parameter_table$se)))
  }
})

test_that("dynamic robust covariance uses free coordinates then delta method", {
  model <- MCMmodel(
    make_moment_summary(), kernel = "dynamic", gaussian_residual = FALSE
  )
  model <- set_dynamic_parameters(model, c(
    phi_A = 0.45, A_lag_to_B = -0.12,
    B_lag_to_A = 0.18, phi_B = 0.35
  ))
  for (variable in c("A", "B")) {
    shape <- paste0("shape_", variable)
    sign <- paste0("sign_", variable)
    model <- MCMparameter(
      model, shape, "free", start = if (variable == "A") 4 else 5,
      transform = "positive"
    )
    model <- MCMparameter(
      model, sign, "fixed", value = if (variable == "A") -1 else 1
    )
  }
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ sign_A * 2 / sqrt(shape_A))
  model <- MCMparameter(model, "kappa_A", "derived",
                        expression = ~ 6 / shape_A)
  model <- MCMparameter(model, "tau_B", "derived",
                        expression = ~ sign_B * 2 / sqrt(shape_B))
  model <- MCMparameter(model, "kappa_B", "derived",
                        expression = ~ 6 / shape_B)
  population <- prepare_population_vcov(
    make_population_summary(model), diag(12) / 10000
  )
  weight <- MCMSEM:::.dynamic_weight_specification(
    population, "diagonal", weight_ridge = 0, require_vcov = TRUE
  )
  covariance <- MCMSEM:::.dynamic_asymptotic_se(
    model, population, weight, se_correction = "robust", weight_ridge = 0
  )
  expect_equal(covariance$jacobian_rank, length(model$param_names))
  expect_identical(dim(covariance$vcov_optimizer),
                   c(length(model$param_names), length(model$param_names)))
  expect_identical(dim(covariance$vcov),
                   c(nrow(model$parameter_table), nrow(model$parameter_table)))
  expect_equal(covariance$se[["sign_A"]], 0)
  expect_true(all(is.finite(covariance$se)))
})

test_that("derived parameters work in dynamic B, cumulant, and Cholesky cells", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "driver", "free", start = 0.2,
                        transform = "bounded", lower = 0, upper = 0.8)
  model <- MCMparameter(model, "phi_A", "derived", expression = ~ driver / 2)
  model <- MCMparameter(model, "tau_A", "derived", expression = ~ -driver)
  model <- MCMparameter(model, "kappa_A", "derived", expression = ~ driver^2)
  model <- MCMparameter(model, "log_sd_G_A", "derived",
                        expression = ~ log(driver))
  implied <- MCMimpliedmoments(model)
  expect_equal(implied$B[1, 1], 0.1)
  expect_equal(model$num_matrices$Tau[1, 1], -0.2)
  expect_equal(model$num_matrices$Kappa[1, 1], 0.04)
  expect_equal(model$num_matrices$L_G[1, 1], log(0.2))
  expect_true(all(is.finite(implied$M4)))
})

test_that("constrained results report all types and can be refitted", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "shape_A", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign_A", "fixed", value = -1)
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ sign_A * 2 / sqrt(shape_A))
  model <- MCMparameter(model, "kappa_A", "derived",
                        expression = ~ 6 / shape_A)
  population <- make_population_summary(model)
  fit <- MCMfit(
    model, population, compute_se = FALSE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1
  )
  expect_named(fit$df, fit$model$parameter_table$name)
  expect_equal(fit$parameter_table$type[fit$parameter_table$parameter == "sign_A"],
               "fixed")
  refit <- MCMfit(
    fit, population, compute_se = FALSE,
    optimizers = "lbfgs", optim_iters = 1, learning_rate = 0.1
  )
  expect_lt(refit$loss, 1e-20)
  expect_equal(MCMparameters(refit), MCMparameters(fit), tolerance = 1e-8)
})

test_that("identification diagnostics are with respect to optimizer coordinates", {
  model <- MCMmodel(make_moment_summary(), kernel = "dynamic")
  model <- MCMparameter(model, "shape_A", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "tau_A", "derived",
                        expression = ~ -2 / sqrt(shape_A))
  model <- MCMparameter(model, "kappa_A", "derived",
                        expression = ~ 6 / shape_A)
  diagnostics <- MCMdiagnostics(model, jacobian = TRUE)
  expect_equal(diagnostics$observed_moment_count, 12)
  expect_equal(diagnostics$independent_free_parameter_count,
               length(model$param_names))
  expect_equal(diagnostics$nominal_df,
               12 - length(model$param_names))
  expect_named(diagnostics$derived_constraints, c("name", "expression"))
  expect_identical(colnames(diagnostics$jacobian), model$param_names)
  expect_true(is.finite(diagnostics$jacobian_condition) ||
                is.infinite(diagnostics$jacobian_condition))
  expect_equal(rownames(diagnostics$near_null_directions), model$param_names)
})

test_that("contemporaneous constrained results can be refitted", {
  model <- MCMmodel(make_moment_summary())
  model <- MCMparameter(model, "shape", "free", start = 4,
                        transform = "positive")
  model <- MCMparameter(model, "sign", "fixed", value = -1)
  model <- MCMparameter(model, "sk1", "derived",
                        expression = ~ sign * 2 / sqrt(shape))
  model <- MCMparameter(model, "k1", "derived",
                        expression = ~ 3 + 6 / shape)
  fit <- MCMfit(
    model, make_moment_summary(), compute_se = FALSE,
    optimizers = "rprop", optim_iters = 1, learning_rate = 0.01
  )
  refit <- MCMfit(
    fit, make_moment_summary(), compute_se = FALSE,
    optimizers = "rprop", optim_iters = 1, learning_rate = 0.01
  )
  expect_identical(refit$model$parameter_table$type,
                   fit$model$parameter_table$type)
  expect_named(refit$df, refit$model$parameter_table$name)
  expect_silent(summary(refit))
})

test_that("legacy graph migration preserves syntax, ordering, and aliases", {
  data <- make_moment_summary()
  default <- MCMmodel(data)
  explicit <- MCMmodel(data, kernel = "contemporaneous")
  alias <- MCMmodel(data, kernel = "static")
  expect_identical(default$param_names, explicit$param_names)
  expect_identical(default$param_names, alias$param_names)
  expect_equal(MCMimpliedmoments(default)$M2, MCMimpliedmoments(explicit)$M2,
               tolerance = 0)
  expect_true(all(MCMparameters(default)$type == "free"))
  old_style <- default$copy()
  old_style$meta_data$kernel <- NULL
  expect_identical(MCMSEM:::.model_kernel(old_style), "contemporaneous")
})
