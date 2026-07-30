.dynamic_diagonal_selector <- function(p, order) {
  out <- matrix(0, p^order, p)
  powers <- p^((order - 1L):0L)
  for (i in seq_len(p)) {
    idx <- 1L + (i - 1L) * sum(powers)
    out[idx, i] <- 1
  }
  out
}

.dynamic_fourth_pairing_selectors <- function(p) {
  grid <- expand.grid(
    l = seq_len(p), k = seq_len(p), j = seq_len(p), i = seq_len(p)
  )
  grid <- as.matrix(grid[, c("i", "j", "k", "l")])
  selector <- function(a, b) {
    out <- matrix(0, nrow(grid), p^2)
    idx <- (grid[, a] - 1L) * p + grid[, b]
    out[cbind(seq_len(nrow(grid)), idx)] <- 1
    out
  }
  list(
    ij = selector(1L, 2L), kl = selector(3L, 4L),
    ik = selector(1L, 3L), jl = selector(2L, 4L),
    il = selector(1L, 4L), jk = selector(2L, 3L)
  )
}

.get_dynamic_torch_components <- function(model, device, dtype) {
  if (.parameter_graph_active(model)) {
    coordinates <- .parameter_optimizer_coordinates(model)
    coordinate_tensor <- torch_tensor(
      unname(coordinates), requires_grad = length(coordinates) > 0L,
      device = device, dtype = dtype
    )
    optimizer_bounds <- .parameter_optimizer_bounds(model)
    p <- model$meta_data$n_phenotypes
    selectors <- list(
      D2 = torch_tensor(.dynamic_diagonal_selector(p, 2L), device = device, dtype = dtype),
      D3 = torch_tensor(.dynamic_diagonal_selector(p, 3L), device = device, dtype = dtype),
      D4 = torch_tensor(.dynamic_diagonal_selector(p, 4L), device = device, dtype = dtype)
    )
    selectors$pairing <- lapply(
      .dynamic_fourth_pairing_selectors(p), torch_tensor,
      device = device, dtype = dtype
    )
    return(list(
      param_names = list(graph = model$param_names),
      par_list = list(graph = coordinate_tensor),
      maps = list(), base_matrices = list(),
      lower = list(graph = torch_tensor(
        unname(optimizer_bounds$L), device = device, dtype = dtype
      )),
      upper = list(graph = torch_tensor(
        unname(optimizer_bounds$U), device = device, dtype = dtype
      )),
      selectors = selectors, p = p,
      gaussian_residual = isTRUE(model$meta_data$gaussian_residual),
      residual_family = .dynamic_residual_family(model),
      dtype = dtype, device = device, parameter_graph = model$copy()
    ))
  }
  matrix_names <- c("B", "Tau", "Kappa", "L_G")
  dtypes <- stats::setNames(rep(list(dtype), length(matrix_names)), matrix_names)
  torch_matrices <- lapply(matrix_names, function(nm) {
    torch_tensor(model$num_matrices[[nm]], device = device, dtype = dtypes[[nm]])
  })
  names(torch_matrices) <- matrix_names
  maps <- stats::setNames(lapply(matrix_names, function(x) list()), matrix_names)
  masks <- lapply(torch_matrices, function(x) torch_ones_like(x, device = device, dtype = dtype))
  param_list <- stats::setNames(vector("list", length(matrix_names)), matrix_names)
  par_values <- stats::setNames(vector("list", length(matrix_names)), matrix_names)

  for (i in seq_along(model$param_coords)) {
    coord <- model$param_coords[[i]]
    matrix_name <- coord[[1]]
    if (!(matrix_name %in% matrix_names)) next
    linear_coords <- coord[[2]]
    multipliers <- coord[[3]]
    for (j in seq_along(linear_coords)) {
      rc <- .r_1to2d_idx(linear_coords[j], nrow(model$num_matrices[[matrix_name]]))
      new_map <- torch_zeros_like(torch_matrices[[matrix_name]], device = device, dtype = dtype)
      new_map[rc[1], rc[2]] <- multipliers[j]
      masks[[matrix_name]][rc[1], rc[2]] <- 0
      parameter_name <- gsub("-", "", model$named_matrices[[matrix_name]][rc[1], rc[2]])
      if (parameter_name %in% names(maps[[matrix_name]])) {
        maps[[matrix_name]][[parameter_name]] <- maps[[matrix_name]][[parameter_name]] + new_map
      } else {
        maps[[matrix_name]][[parameter_name]] <- new_map
        param_list[[matrix_name]] <- c(param_list[[matrix_name]], parameter_name)
        par_values[[matrix_name]] <- c(
          par_values[[matrix_name]],
          model$param_values[model$param_names == parameter_name][1]
        )
      }
    }
  }

  for (nm in matrix_names) {
    maps[[nm]] <- if (length(maps[[nm]]) > 0L) {
      torch_dstack(maps[[nm]])
    } else {
      torch_zeros_like(torch_matrices[[nm]], device = device, dtype = dtype)
    }
    if (length(maps[[nm]]$shape) == 2L) {
      shape <- as.numeric(maps[[nm]]$shape)
      maps[[nm]] <- torch_reshape(maps[[nm]], c(shape, 1L))
    }
  }

  base_matrices <- lapply(matrix_names, function(nm) {
    torch_matrices[[nm]] * masks[[nm]]
  })
  names(base_matrices) <- matrix_names
  par_list <- lapply(matrix_names, function(nm) {
    if (length(param_list[[nm]]) > 0L) {
      torch_tensor(
        as.numeric(par_values[[nm]]), requires_grad = TRUE,
        device = device, dtype = dtype
      )
    } else {
      torch_tensor(1, device = device, dtype = dtype)
    }
  })
  names(par_list) <- matrix_names

  lower <- upper <- list()
  for (nm in matrix_names) {
    if (length(param_list[[nm]]) > 0L) {
      lower[[nm]] <- torch_tensor(
        as.numeric(model$bounds["L", param_list[[nm]], drop = TRUE]),
        device = device, dtype = dtype
      )
      upper[[nm]] <- torch_tensor(
        as.numeric(model$bounds["U", param_list[[nm]], drop = TRUE]),
        device = device, dtype = dtype
      )
    }
  }

  p <- model$meta_data$n_phenotypes
  selectors <- list(
    D2 = torch_tensor(.dynamic_diagonal_selector(p, 2L), device = device, dtype = dtype),
    D3 = torch_tensor(.dynamic_diagonal_selector(p, 3L), device = device, dtype = dtype),
    D4 = torch_tensor(.dynamic_diagonal_selector(p, 4L), device = device, dtype = dtype)
  )
  pairing <- .dynamic_fourth_pairing_selectors(p)
  selectors$pairing <- lapply(pairing, torch_tensor, device = device, dtype = dtype)

  list(
    param_names = param_list,
    par_list = par_list,
    maps = maps,
    base_matrices = base_matrices,
    lower = lower,
    upper = upper,
    selectors = selectors,
    p = p,
    gaussian_residual = isTRUE(model$meta_data$gaussian_residual),
    residual_family = .dynamic_residual_family(model),
    dtype = dtype,
    device = device
  )
}

.dynamic_torch_matrix <- function(components, name) {
  if (!is.null(components$parameter_graph)) {
    values <- .parameter_values_torch(
      components$parameter_graph, components$par_list$graph
    )
    return(.parameter_torch_matrix(
      components$parameter_graph, values, name,
      components$par_list$graph
    ))
  }
  components$base_matrices[[name]] + torch_sum(
    components$maps[[name]] * components$par_list[[name]], dim = 3
  )
}

.dynamic_torch_kron_power <- function(B, order) {
  out <- B
  if (order > 1L) {
    for (i in 2:order) out <- .torch_kron(out, B)
  }
  out
}

.dynamic_torch_cumulant <- function(B, innovation_values, selector,
                                      order, p) {
  transition <- .dynamic_torch_kron_power(B, order)
  lhs <- torch_eye(p^order, device = B$device, dtype = B$dtype) - transition
  D <- torch_matmul(selector, innovation_values)
  linalg_solve(lhs, D)
}

.dynamic_torch_gaussian_covariance <- function(raw_L) {
  L <- torch_tril(raw_L, diagonal = -1L) + torch_diag(torch_exp(torch_diag(raw_L)))
  list(L = L, Psi_G = torch_matmul(L, torch_transpose(L, 1, 2)))
}

.dynamic_torch_vector_power <- function(x, order) {
  out <- x
  if (order > 1L) {
    for (i in 2:order) out <- .torch_kron(out, x)
  }
  out
}

.dynamic_torch_residual_cumulants <- function(components, reference) {
  p <- components$p
  family <- components$residual_family
  zero2 <- torch_zeros(c(p, p), device = reference$device,
                       dtype = reference$dtype)
  zero3 <- torch_zeros(p^3, device = reference$device,
                       dtype = reference$dtype)
  zero4 <- torch_zeros(p^4, device = reference$device,
                       dtype = reference$dtype)
  if (identical(family, "none")) {
    return(list(
      M2 = zero2, C3 = zero3, K4 = zero4, Psi_G = zero2, L_G = zero2,
      loadings = torch_zeros(p, device = reference$device,
                             dtype = reference$dtype),
      shape = torch_tensor(NaN, device = reference$device,
                           dtype = reference$dtype),
      skewness = torch_tensor(0, device = reference$device,
                              dtype = reference$dtype),
      excess_kurtosis = torch_tensor(0, device = reference$device,
                                     dtype = reference$dtype)
    ))
  }
  if (identical(family, "gaussian")) {
    gaussian <- .dynamic_torch_gaussian_covariance(
      .dynamic_torch_matrix(components, "L_G")
    )
    return(list(
      M2 = gaussian$Psi_G, C3 = zero3, K4 = zero4,
      Psi_G = gaussian$Psi_G, L_G = gaussian$L,
      loadings = torch_zeros(p, device = reference$device,
                             dtype = reference$dtype),
      shape = torch_tensor(NaN, device = reference$device,
                           dtype = reference$dtype),
      skewness = torch_tensor(0, device = reference$device,
                              dtype = reference$dtype),
      excess_kurtosis = torch_tensor(0, device = reference$device,
                                     dtype = reference$dtype)
    ))
  }
  if (!identical(family, "common_gamma")) {
    stop("Unsupported dynamic residual family `", family, "`.", call. = FALSE)
  }
  loadings <- torch_flatten(
    .dynamic_torch_matrix(components, "Lambda_Gamma")
  )
  shape <- torch_flatten(
    .dynamic_torch_matrix(components, "Shape_Gamma")
  )[1]
  skewness <- 2 / torch_sqrt(shape)
  excess_kurtosis <- 6 / shape
  loading_column <- torch_reshape(loadings, c(p, 1L))
  M2 <- torch_matmul(loading_column, torch_transpose(loading_column, 1L, 2L))
  list(
    M2 = M2,
    C3 = skewness * .dynamic_torch_vector_power(loadings, 3L),
    K4 = excess_kurtosis * .dynamic_torch_vector_power(loadings, 4L),
    Psi_G = zero2, L_G = zero2, loadings = loadings, shape = shape,
    skewness = skewness, excess_kurtosis = excess_kurtosis
  )
}

.dynamic_torch_raw_fourth <- function(K4_vector, Sigma, selectors, p) {
  sigma <- torch_flatten(Sigma)
  select <- function(name) torch_matmul(selectors[[name]], sigma)
  pairings <- select("ij") * select("kl") +
    select("ik") * select("jl") +
    select("il") * select("jk")
  torch_reshape(K4_vector + pairings, c(p, p^3))
}

.get_dynamic_predicted_matrices <- function(components) {
  p <- components$p
  B <- .dynamic_torch_matrix(components, "B")
  tau <- torch_flatten(.dynamic_torch_matrix(components, "Tau"))
  kappa <- torch_flatten(.dynamic_torch_matrix(components, "Kappa"))
  ones <- torch_ones(p, device = B$device, dtype = B$dtype)

  C2_vector <- .dynamic_torch_cumulant(
    B, ones, components$selectors$D2, 2L, p
  )
  C3_vector <- .dynamic_torch_cumulant(
    B, tau, components$selectors$D3, 3L, p
  )
  K4_vector <- .dynamic_torch_cumulant(
    B, kappa, components$selectors$D4, 4L, p
  )

  within_M2 <- torch_reshape(C2_vector, c(p, p))
  residual <- .dynamic_torch_residual_cumulants(components, B)
  total_C3_vector <- C3_vector + residual$C3
  total_K4_vector <- K4_vector + residual$K4
  Sigma <- within_M2 + residual$M2
  M3 <- torch_reshape(total_C3_vector, c(p, p^2))
  K4 <- torch_reshape(total_K4_vector, c(p, p^3))
  M4 <- .dynamic_torch_raw_fourth(
    total_K4_vector, Sigma, components$selectors$pairing, p
  )
  rho <- torch_max(torch_abs(linalg_eigvals(B)))

  list(
    M2 = Sigma, M3 = M3, M4 = M4, K4 = K4,
    within_M2 = within_M2,
    dynamic_M3 = torch_reshape(C3_vector, c(p, p^2)),
    dynamic_K4 = torch_reshape(K4_vector, c(p, p^3)),
    residual_M2 = residual$M2,
    residual_M3 = torch_reshape(residual$C3, c(p, p^2)),
    residual_K4 = torch_reshape(residual$K4, c(p, p^3)),
    Psi_G = residual$Psi_G, L_G = residual$L_G,
    residual_loadings = residual$loadings,
    residual_shape = residual$shape,
    residual_skewness = residual$skewness,
    residual_excess_kurtosis = residual$excess_kurtosis,
    B = B, tau = tau, kappa = kappa, spectral_radius = rho
  )
}

.dynamic_grad_parameters <- function(components) {
  components$par_list[vapply(
    components$par_list,
    function(x) isTRUE(x$requires_grad),
    logical(1)
  )]
}

.dynamic_bounds_penalty <- function(components) {
  parameters <- .dynamic_grad_parameters(components)
  if (length(parameters) == 0L) {
    return(torch_tensor(0, device = components$device, dtype = components$dtype))
  }
  par <- torch_cat(parameters)
  lower <- torch_cat(components$lower[names(parameters)])
  upper <- torch_cat(components$upper[names(parameters)])
  torch_sum(torch_square(torch_relu(lower - par))) +
    torch_sum(torch_square(torch_relu(par - upper)))
}

.dynamic_objective <- function(components, lossfunc, m2v_masks,
                               observed, use_bounds,
                               outofbounds_penalty,
                               stationarity_limit,
                               stationarity_penalty,
                               moment_weight, moment_grids,
                               loss_type) {
  predicted <- .get_dynamic_predicted_matrices(components)
  moment_loss <- if (identical(loss_type, "mse")) {
    residual <- .dynamic_torch_moment_vector(predicted, moment_grids) -
      .dynamic_torch_moment_vector(observed, moment_grids)
    torch_dot(residual, torch_matmul(moment_weight, residual))
  } else {
    .calc_loss(
      lossfunc, predicted, m2v_masks,
      observed$M2, observed$M3, observed$M4,
      use_skewness = TRUE, use_kurtosis = TRUE
    )
  }
  penalty <- stationarity_penalty * torch_square(
    torch_relu(predicted$spectral_radius - stationarity_limit)
  )
  if (isTRUE(use_bounds)) {
    penalty <- penalty + outofbounds_penalty * .dynamic_bounds_penalty(components)
  }
  moment_loss + penalty
}

.dynamic_parameter_vector <- function(components) {
  if (!is.null(components$parameter_graph)) {
    coordinates <- as.numeric(torch_tensor(
      components$par_list$graph, device = torch_device("cpu")
    ))
    values <- .parameter_values_base(
      components$parameter_graph, coordinates, optimizer_scale = TRUE
    )
    return(unname(values[components$parameter_graph$param_names]))
  }
  parameters <- .dynamic_grad_parameters(components)
  as.numeric(torch_tensor(
    torch_cat(parameters), device = torch_device("cpu")
  ))
}

.dynamic_gradient_snapshot <- function(components) {
  out <- list()
  for (nm in names(components$par_list)) {
    par <- components$par_list[[nm]]
    if (isTRUE(par$requires_grad) && !is.null(par$grad)) {
      values <- as.numeric(torch_tensor(par$grad, device = torch_device("cpu")))
      names(values) <- components$param_names[[nm]]
      out[[nm]] <- values
    }
  }
  out
}

.torch_fit_dynamic <- function(model, observed, m2v_masks, optimizers,
                               optim_iters, learning_rate, lossfunc,
                               device, dtype, use_bounds,
                               outofbounds_penalty, stationarity_limit,
                               stationarity_penalty, moment_weight,
                               moment_grids, loss_type,
                               monitor_grads = FALSE,
                               verbose = FALSE) {
  components <- .get_dynamic_torch_components(model, device, dtype)
  parameters <- .dynamic_grad_parameters(components)
  loss_history <- numeric()
  gradient_history <- list()

  objective <- function() {
    .dynamic_objective(
      components, lossfunc, m2v_masks, observed, use_bounds,
      outofbounds_penalty, stationarity_limit, stationarity_penalty,
      moment_weight, moment_grids, loss_type
    )
  }

  initial_loss <- as.numeric(torch_tensor(objective(), device = torch_device("cpu")))
  if (length(parameters) > 0L) for (noptim in seq_along(optimizers)) {
    optimizer_name <- optimizers[noptim]
    if (isTRUE(verbose)) {
      cat(sprintf("  optimizer=%s, learning_rate=%g\n",
                  optimizer_name, learning_rate[noptim]))
    }
    if (!identical(optimizer_name, "lbfgs")) {
      optimizer <- .get_optimfunc(optimizer_name)(
        parameters, lr = learning_rate[noptim]
      )
      for (iteration in seq_len(optim_iters[noptim])) {
        optimizer$zero_grad()
        loss <- objective()
        loss$backward()
        if (isTRUE(monitor_grads)) {
          gradient_history[[length(gradient_history) + 1L]] <-
            .dynamic_gradient_snapshot(components)
        }
        loss_history <- c(
          loss_history,
          as.numeric(torch_tensor(loss$detach(), device = torch_device("cpu")))
        )
        optimizer$step()
        if (isTRUE(verbose)) {
          cat(sprintf("\r    iteration=%d, penalized loss=%.8g", iteration,
                      tail(loss_history, 1L)))
          flush.console()
        }
      }
      if (isTRUE(verbose)) cat("\n")
    } else {
      optimizer <- optim_lbfgs(parameters, lr = learning_rate[noptim])
      closure <- function() {
        optimizer$zero_grad()
        loss <- objective()
        loss$backward()
        if (isTRUE(monitor_grads)) {
          gradient_history[[length(gradient_history) + 1L]] <<-
            .dynamic_gradient_snapshot(components)
        }
        loss_history <<- c(
          loss_history,
          as.numeric(torch_tensor(loss$detach(), device = torch_device("cpu")))
        )
        loss
      }
      for (iteration in seq_len(optim_iters[noptim])) {
        optimizer$step(closure)
        if (isTRUE(verbose)) {
          cat(sprintf("\r    iteration=%d, penalized loss=%.8g", iteration,
                      tail(loss_history, 1L)))
          flush.console()
        }
      }
      if (isTRUE(verbose)) cat("\n")
    }
  }

  predicted <- .get_dynamic_predicted_matrices(components)
  final_loss <- if (identical(loss_type, "mse")) {
    residual <- .dynamic_torch_moment_vector(predicted, moment_grids) -
      .dynamic_torch_moment_vector(observed, moment_grids)
    torch_dot(residual, torch_matmul(moment_weight, residual))
  } else {
    .calc_loss(
      lossfunc, predicted, m2v_masks,
      observed$M2, observed$M3, observed$M4, TRUE, TRUE
    )
  }
  if (!isTRUE(monitor_grads) && length(parameters) > 0L) {
    for (parameter in parameters) {
      if (!is.null(parameter$grad)) parameter$grad$zero_()
    }
    objective()$backward()
    gradient_history[[1L]] <- .dynamic_gradient_snapshot(components)
  }
  list(
    parameters = .dynamic_parameter_vector(components),
    optimizer_coordinates = if (!is.null(components$parameter_graph)) {
      as.numeric(torch_tensor(
        components$par_list$graph, device = torch_device("cpu")
      ))
    } else {
      .dynamic_parameter_vector(components)
    },
    predicted = predicted,
    initial_loss = initial_loss,
    loss = as.numeric(torch_tensor(final_loss, device = torch_device("cpu"))),
    loss_history = loss_history,
    gradient_history = gradient_history,
    last_gradient = .dynamic_gradient_snapshot(components)
  )
}
