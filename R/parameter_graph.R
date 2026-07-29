# Canonical free, fixed, and derived parameter graph ----------------------

.empty_parameter_table <- function() {
  out <- data.frame(
    name = character(), type = character(), value = numeric(),
    start = numeric(), lower = numeric(), upper = numeric(),
    transform = character(), expression = character(),
    optimizer_index = integer(), optimizer_value = numeric(),
    auxiliary = logical(), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$dependencies <- I(list())
  out$matrix_locations <- I(list())
  out
}

.parameter_table_copy <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(.empty_parameter_table())
  out <- x
  if ("dependencies" %in% names(out)) {
    out$dependencies <- I(lapply(out$dependencies, function(z) as.character(z)))
  }
  if ("matrix_locations" %in% names(out)) {
    out$matrix_locations <- I(lapply(out$matrix_locations, function(z) {
      if (is.data.frame(z)) z[, , drop = FALSE] else z
    }))
  }
  out
}

.parameter_table_row <- function(name, type = "free", value = NA_real_,
                                 start = NA_real_, lower = -Inf,
                                 upper = Inf, transform = "identity",
                                 expression = "", dependencies = character(),
                                 optimizer_index = NA_integer_,
                                 optimizer_value = NA_real_, auxiliary = FALSE,
                                 matrix_locations = data.frame()) {
  out <- data.frame(
    name = as.character(name), type = as.character(type),
    value = as.numeric(value), start = as.numeric(start),
    lower = as.numeric(lower), upper = as.numeric(upper),
    transform = as.character(transform), expression = as.character(expression),
    optimizer_index = as.integer(optimizer_index),
    optimizer_value = as.numeric(optimizer_value),
    auxiliary = isTRUE(auxiliary), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$dependencies <- I(list(as.character(dependencies)))
  out$matrix_locations <- I(list(matrix_locations))
  out
}

.parameter_field <- function(model, name, default) {
  tryCatch(model[[name]], error = function(e) default)
}

.parameter_label_parts <- function(x) {
  negative <- startsWith(x, "-")
  list(name = ifelse(negative, substring(x, 2L), x),
       multiplier = ifelse(negative, -1, 1))
}

.parameter_matrix_locations <- function(model) {
  order <- character()
  locations <- list()
  for (matrix_name in names(model$named_matrices)) {
    labels <- as.vector(model$named_matrices[[matrix_name]])
    numeric_label <- !is.na(suppressWarnings(as.numeric(labels)))
    for (index in which(!numeric_label)) {
      parts <- .parameter_label_parts(labels[index])
      name <- parts$name
      if (!nzchar(name)) next
      if (!(name %in% order)) order <- c(order, name)
      rc <- .r_1to2d_idx(index, nrow(model$named_matrices[[matrix_name]]))
      location <- data.frame(
        matrix = matrix_name, index = index, row = rc[1], col = rc[2],
        multiplier = parts$multiplier, stringsAsFactors = FALSE
      )
      locations[[name]] <- rbind(locations[[name]], location)
    }
  }
  list(order = order, locations = locations)
}

.parameter_default_bounds <- function(model, name, location) {
  existing <- model$bounds
  if (is.data.frame(existing) && name %in% colnames(existing)) {
    value <- as.numeric(existing[c("L", "U"), name, drop = TRUE])
    names(value) <- c("L", "U")
    if (!all(is.na(value))) return(value)
  }
  if (is.data.frame(location) && nrow(location)) {
    matrix_name <- location$matrix[1L]
    if (identical(.model_kernel(model), "dynamic")) {
      if (matrix_name == "B") {
        diagonal <- all(location$row == location$col)
        return(c(L = if (diagonal) 0 else -0.98, U = 0.98))
      }
      if (matrix_name == "Tau") return(c(L = -100, U = 100))
      if (matrix_name == "Kappa") return(c(L = -1.99, U = 100))
      if (matrix_name == "L_G") {
        diagonal <- all(location$row == location$col)
        return(c(L = if (diagonal) log(1e-6) else -100,
                 U = if (diagonal) log(1e2) else 100))
      }
    } else {
      defaults <- model$meta_data$bound_default
      if (is.null(defaults)) defaults <- model$meta_data$bound_defaults
      prefix <- sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", name))
      if (!is.null(defaults$L[[prefix]]) && !is.null(defaults$U[[prefix]])) {
        return(c(L = defaults$L[[prefix]], U = defaults$U[[prefix]]))
      }
    }
  }
  value <- tryCatch(.model_default_bounds(model, name), error = function(e) NULL)
  if (is.null(value) || length(value) != 2L || all(is.na(value))) {
    value <- c(L = -Inf, U = Inf)
  }
  value[is.na(value)] <- c(-Inf, Inf)[is.na(value)]
  value
}

.parameter_expression_call <- function(expression) {
  if (inherits(expression, "formula")) {
    if (length(expression) != 2L) {
      stop("A derived parameter expression must be a one-sided formula.", call. = FALSE)
    }
    return(expression[[2L]])
  }
  if (is.call(expression) || is.name(expression) ||
      (is.numeric(expression) && length(expression) == 1L)) {
    return(expression)
  }
  stop("`expression` must be a one-sided formula or a validated language expression.",
       call. = FALSE)
}

.parameter_expression_dependencies <- function(node) {
  if (is.numeric(node)) {
    if (length(node) != 1L || !is.finite(node)) {
      stop("Expression constants must be finite scalar numbers.", call. = FALSE)
    }
    return(character())
  }
  if (is.name(node)) return(as.character(node))
  if (!is.call(node)) {
    stop("Unsupported object in a derived parameter expression.", call. = FALSE)
  }
  fun <- as.character(node[[1L]])
  binary <- c("+", "-", "*", "/", "^")
  unary_functions <- c("sqrt", "exp", "log", "softplus", "logistic")
  nargs <- length(node) - 1L
  if (fun %in% binary) {
    valid_arity <- if (fun %in% c("+", "-")) nargs %in% c(1L, 2L) else nargs == 2L
    if (!valid_arity) {
      stop("Operator `", fun, "` has an invalid number of arguments.", call. = FALSE)
    }
  } else if (fun %in% unary_functions) {
    if (nargs != 1L) {
      stop("Function `", fun, "` requires exactly one argument.", call. = FALSE)
    }
  } else {
    stop("Unsupported function or operator `", fun, "` in derived expression.",
         call. = FALSE)
  }
  unique(unlist(lapply(as.list(node)[-1L], .parameter_expression_dependencies),
                use.names = FALSE))
}

.parameter_topological_order <- function(table) {
  if (nrow(table) == 0L) return(character())
  derived <- table$name[table$type == "derived"]
  state <- stats::setNames(integer(length(derived)), derived)
  out <- character()
  visit <- function(name, path = character()) {
    if (state[[name]] == 2L) return(invisible(NULL))
    if (state[[name]] == 1L) {
      cycle <- c(path[match(name, path):length(path)], name)
      stop("Derived parameter dependency cycle: ", paste(cycle, collapse = " -> "),
           call. = FALSE)
    }
    state[[name]] <<- 1L
    deps <- table$dependencies[[match(name, table$name)]]
    for (dependency in intersect(deps, derived)) visit(dependency, c(path, name))
    state[[name]] <<- 2L
    out <<- c(out, name)
    invisible(NULL)
  }
  for (name in derived) visit(name)
  out
}

.parameter_validate_table <- function(table, expressions) {
  if (anyDuplicated(table$name)) {
    stop("Parameter names must be unique.", call. = FALSE)
  }
  if (any(!nzchar(table$name)) || anyNA(table$name)) {
    stop("Parameter names must be non-empty and non-missing.", call. = FALSE)
  }
  if (any(!(table$type %in% c("free", "fixed", "derived")))) {
    stop("Parameter type must be `free`, `fixed`, or `derived`.", call. = FALSE)
  }
  for (i in seq_len(nrow(table))) {
    name <- table$name[i]
    type <- table$type[i]
    if (type == "free") {
      if (!is.finite(table$value[i]) || !is.finite(table$start[i])) {
        stop("Free parameter `", name, "` requires finite current and starting values.",
             call. = FALSE)
      }
      if (!(table$transform[i] %in% c("identity", "positive", "bounded"))) {
        stop("Unknown transformation for free parameter `", name, "`.", call. = FALSE)
      }
      if (table$lower[i] >= table$upper[i]) {
        stop("Lower bound must be below upper bound for parameter `", name, "`.",
             call. = FALSE)
      }
      if (table$transform[i] == "positive" &&
          table$value[i] <= table$lower[i]) {
        stop("Positive parameter `", name, "` must start above its lower bound.",
             call. = FALSE)
      }
      if (table$transform[i] == "bounded" &&
          (!is.finite(table$lower[i]) || !is.finite(table$upper[i]) ||
           table$value[i] <= table$lower[i] || table$value[i] >= table$upper[i])) {
        stop("Bounded parameter `", name,
             "` must have finite bounds and start strictly between them.", call. = FALSE)
      }
    } else if (type == "fixed") {
      if (!is.finite(table$value[i])) {
        stop("Fixed parameter `", name, "` must have a finite value.", call. = FALSE)
      }
    } else {
      if (is.null(expressions[[name]])) {
        stop("Derived parameter `", name, "` has no expression.", call. = FALSE)
      }
      unknown <- setdiff(table$dependencies[[i]], table$name)
      if (length(unknown)) {
        stop("Unknown symbol", if (length(unknown) > 1L) "s" else "", " in `", name,
             "`: ", paste(unknown, collapse = ", "), call. = FALSE)
      }
    }
  }
  .parameter_topological_order(table)
}

.parameter_softplus <- function(x) log1p(exp(-abs(x))) + pmax(x, 0)

.parameter_inverse_softplus <- function(x) {
  ifelse(x > 20, x + log1p(-exp(-x)), log(expm1(x)))
}

.parameter_natural_to_optimizer_value <- function(value, transform, lower, upper) {
  switch(
    transform,
    identity = value,
    positive = .parameter_inverse_softplus(value - lower),
    bounded = stats::qlogis((value - lower) / (upper - lower)),
    stop("Unsupported parameter transformation.", call. = FALSE)
  )
}

.parameter_optimizer_to_natural_value <- function(value, transform, lower, upper) {
  switch(
    transform,
    identity = value,
    positive = lower + .parameter_softplus(value),
    bounded = lower + (upper - lower) * stats::plogis(value),
    stop("Unsupported parameter transformation.", call. = FALSE)
  )
}

.parameter_eval_base_node <- function(node, values) {
  if (is.numeric(node)) return(as.numeric(node))
  if (is.name(node)) return(values[[as.character(node)]])
  fun <- as.character(node[[1L]])
  args <- lapply(as.list(node)[-1L], .parameter_eval_base_node, values = values)
  if (fun == "+") return(if (length(args) == 1L) args[[1L]] else args[[1L]] + args[[2L]])
  if (fun == "-") return(if (length(args) == 1L) -args[[1L]] else args[[1L]] - args[[2L]])
  if (fun == "*") return(args[[1L]] * args[[2L]])
  if (fun == "/") return(args[[1L]] / args[[2L]])
  if (fun == "^") return(args[[1L]] ^ args[[2L]])
  if (fun == "sqrt") return(sqrt(args[[1L]]))
  if (fun == "exp") return(exp(args[[1L]]))
  if (fun == "log") return(log(args[[1L]]))
  if (fun == "softplus") return(.parameter_softplus(args[[1L]]))
  if (fun == "logistic") return(stats::plogis(args[[1L]]))
  stop("Internal error: unvalidated expression node.", call. = FALSE)
}

.parameter_values_base <- function(model, parameters = NULL,
                                   optimizer_scale = FALSE,
                                   check_finite = TRUE) {
  table <- model$parameter_table
  expressions <- model$parameter_expressions
  free_rows <- which(table$type == "free")
  free_names <- table$name[free_rows]
  if (is.null(parameters)) parameters <- model$param_values
  if (length(parameters) != length(free_rows)) {
    stop("Parameter vector has the wrong length; expected ", length(free_rows),
         " independent free parameters.", call. = FALSE)
  }
  if (!is.null(names(parameters)) && all(free_names %in% names(parameters))) {
    parameters <- parameters[free_names]
  }
  parameters <- as.numeric(parameters)
  values <- list()
  for (j in seq_along(free_rows)) {
    i <- free_rows[j]
    values[[table$name[i]]] <- if (isTRUE(optimizer_scale)) {
      .parameter_optimizer_to_natural_value(
        parameters[j], table$transform[i], table$lower[i], table$upper[i]
      )
    } else parameters[j]
  }
  for (i in which(table$type == "fixed")) values[[table$name[i]]] <- table$value[i]
  for (name in .parameter_topological_order(table)) {
    values[[name]] <- suppressWarnings(
      .parameter_eval_base_node(expressions[[name]], values)
    )
  }
  out <- stats::setNames(
    vapply(table$name, function(name) as.numeric(values[[name]]), numeric(1)),
    table$name
  )
  if (isTRUE(check_finite) && any(!is.finite(out))) {
    bad <- names(out)[!is.finite(out)]
    stop("Parameter graph evaluated to a non-finite value for: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

.parameter_optimizer_coordinates <- function(model, parameters = NULL) {
  table <- model$parameter_table
  free <- which(table$type == "free")
  if (is.null(parameters)) parameters <- model$param_values
  if (!is.null(names(parameters)) && all(table$name[free] %in% names(parameters))) {
    parameters <- parameters[table$name[free]]
  }
  if (length(parameters) != length(free)) {
    stop("Free parameter vector has the wrong length.", call. = FALSE)
  }
  out <- vapply(seq_along(free), function(j) {
    i <- free[j]
    .parameter_natural_to_optimizer_value(
      as.numeric(parameters[j]), table$transform[i], table$lower[i], table$upper[i]
    )
  }, numeric(1))
  names(out) <- table$name[free]
  if (any(!is.finite(out))) {
    stop("Starting values could not be mapped to finite optimizer coordinates.",
         call. = FALSE)
  }
  out
}

.parameter_optimizer_bounds <- function(model, finite = 1e10) {
  table <- model$parameter_table
  free <- which(table$type == "free")
  lower <- upper <- numeric(length(free))
  for (j in seq_along(free)) {
    i <- free[j]
    if (table$transform[i] == "identity") {
      lower[j] <- table$lower[i]
      upper[j] <- table$upper[i]
    } else {
      lower[j] <- -finite
      upper[j] <- finite
    }
  }
  lower[!is.finite(lower)] <- -finite
  upper[!is.finite(upper)] <- finite
  list(L = stats::setNames(lower, table$name[free]),
       U = stats::setNames(upper, table$name[free]))
}

.parameter_graph_apply_to_model <- function(model, parameters = NULL,
                                            optimizer_scale = FALSE) {
  values <- .parameter_values_base(model, parameters, optimizer_scale)
  table <- model$parameter_table
  table$value <- unname(values[table$name])
  free <- table$type == "free"
  if (any(free)) {
    table$optimizer_value[free] <- vapply(which(free), function(i) {
      .parameter_natural_to_optimizer_value(
        table$value[i], table$transform[i], table$lower[i], table$upper[i]
      )
    }, numeric(1))
  }
  if (!isTRUE(optimizer_scale)) model$param_values <- unname(values[table$name[free]])
  model$parameter_table <- table
  for (i in seq_len(nrow(table))) {
    locations <- table$matrix_locations[[i]]
    if (!is.data.frame(locations) || nrow(locations) == 0L) next
    for (j in seq_len(nrow(locations))) {
      matrix_name <- locations$matrix[j]
      model$num_matrices[[matrix_name]][locations$index[j]] <-
        locations$multiplier[j] * values[[table$name[i]]]
    }
  }
  invisible(values)
}

.parameter_graph_rebuild <- function(model) {
  old_table <- .parameter_field(model, "parameter_table", data.frame())
  if (!is.data.frame(old_table) || !all(c("name", "type") %in% names(old_table))) {
    old_table <- .empty_parameter_table()
  }
  old_expressions <- .parameter_field(model, "parameter_expressions", list())
  old_names <- as.character(.parameter_field(model, "param_names", character()))
  old_values <- as.numeric(.parameter_field(model, "param_values", numeric()))
  names(old_values) <- old_names
  old_starts <- tryCatch({
    out <- as.numeric(model$start_values$get("start", model$start_values$getcolnames()))
    names(out) <- model$start_values$getcolnames()
    out
  }, error = function(e) numeric())
  old_bounds <- model$bounds

  scan <- .parameter_matrix_locations(model)
  keep_auxiliary <- if (nrow(old_table)) {
    old_table$name[old_table$auxiliary & !(old_table$name %in% scan$order)]
  } else character()
  target_names <- c(scan$order, setdiff(keep_auxiliary, scan$order))
  table <- .empty_parameter_table()
  expressions <- list()

  for (name in target_names) {
    old_i <- match(name, old_table$name)
    location <- scan$locations[[name]]
    if (is.null(location)) location <- data.frame()
    if (!is.na(old_i)) {
      row <- old_table[old_i, , drop = FALSE]
      row$matrix_locations <- I(list(location))
      if (name %in% names(old_values) && row$type == "free") {
        row$value <- old_values[[name]]
      }
      if (name %in% names(old_starts) && row$type == "free") {
        row$start <- old_starts[[name]]
      }
      if (is.data.frame(old_bounds) && name %in% colnames(old_bounds) &&
          row$type == "free") {
        row$lower <- as.numeric(old_bounds["L", name])
        row$upper <- as.numeric(old_bounds["U", name])
      }
      if (row$type == "derived") expressions[[name]] <- old_expressions[[name]]
    } else {
      if (nrow(location)) {
        first <- location[1L, ]
        initial <- model$num_matrices[[first$matrix]][first$index] /
          first$multiplier
        bounds <- .parameter_default_bounds(model, name, first)
      } else {
        initial <- 0
        bounds <- c(L = -Inf, U = Inf)
      }
      row <- .parameter_table_row(
        name = name, type = "free", value = initial, start = initial,
        lower = bounds[["L"]], upper = bounds[["U"]],
        transform = "identity", auxiliary = !nrow(location),
        matrix_locations = location
      )
    }
    table <- rbind(table, row)
  }

  if (nrow(table)) {
    for (i in which(table$type == "derived")) {
      name <- table$name[i]
      call <- expressions[[name]]
      table$dependencies[[i]] <- .parameter_expression_dependencies(call)
      table$expression[i] <- paste(deparse(call, width.cutoff = 500L), collapse = " ")
    }
    for (i in which(table$type != "derived")) {
      table$dependencies[[i]] <- character()
      table$expression[i] <- ""
    }
  }

  # Build provisional free coordinates before looking up legacy matrix defaults.
  free <- which(table$type == "free")
  model$param_names <- table$name[free]
  model$free_param_names <- table$name[free]
  model$param_values <- table$value[free]
  model$param_coords <- lapply(free, function(i) {
    locations <- table$matrix_locations[[i]]
    if (!is.data.frame(locations) || nrow(locations) == 0L) {
      return(list("", integer(), numeric()))
    }
    first_matrix <- locations$matrix[1L]
    same <- locations$matrix == first_matrix
    list(first_matrix, locations$index[same], locations$multiplier[same])
  })

  for (i in free) {
    if (is.na(table$lower[i]) || is.na(table$upper[i])) {
      defaults <- .parameter_default_bounds(
        model, table$name[i], table$matrix_locations[[i]]
      )
      if (is.na(table$lower[i])) table$lower[i] <- defaults[["L"]]
      if (is.na(table$upper[i])) table$upper[i] <- defaults[["U"]]
    }
  }
  .parameter_validate_table(table, expressions)
  model$parameter_table <- table
  model$parameter_expressions <- expressions

  values <- .parameter_values_base(model, table$value[free])
  table <- model$parameter_table
  table$value <- unname(values[table$name])
  derived_or_fixed <- table$type != "free"
  table$start[derived_or_fixed] <- table$value[derived_or_fixed]
  table$optimizer_index <- NA_integer_
  table$optimizer_value <- NA_real_
  if (length(free)) {
    table$optimizer_index[free] <- seq_along(free)
    table$optimizer_value[free] <- unname(
      .parameter_optimizer_coordinates(model, table$value[free])
    )
  }
  model$parameter_table <- table
  model$param_names <- table$name[free]
  model$free_param_names <- table$name[free]
  model$param_values <- table$value[free]

  starts <- as.data.frame(matrix(table$start[free], nrow = 1L),
                          check.names = FALSE)
  colnames(starts) <- table$name[free]
  rownames(starts) <- "start"
  model$start_values <- mcmstartvaluesclass(starts)
  bounds <- data.frame(row.names = c("L", "U"))
  if (length(free)) {
    bounds <- as.data.frame(rbind(L = table$lower[free], U = table$upper[free]),
                            check.names = FALSE)
    colnames(bounds) <- table$name[free]
  }
  model$bounds <- bounds
  .parameter_graph_apply_to_model(model)
  invisible(model)
}

.ensure_parameter_graph <- function(model) {
  table <- .parameter_field(model, "parameter_table", data.frame())
  if (!is.data.frame(table) || !all(c("name", "type", "matrix_locations") %in% names(table))) {
    model$parse()
  }
  invisible(model)
}

.parameter_graph_active <- function(model) {
  .ensure_parameter_graph(model)
  table <- model$parameter_table
  if (nrow(table) == 0L) return(FALSE)
  if (any(table$type != "free") || any(table$auxiliary) ||
      any(table$transform != "identity") ||
      any(!is.finite(table$lower) | !is.finite(table$upper))) return(TRUE)
  spans <- vapply(table$matrix_locations, function(locations) {
    is.data.frame(locations) && length(unique(locations$matrix)) > 1L
  }, logical(1))
  any(spans)
}

#' Define or convert a model parameter
#'
#' Adds an auxiliary free or fixed parameter, or converts an existing matrix
#' parameter among free, fixed, and safely derived forms. Only free parameters
#' are optimized. Derived expressions support numeric constants, symbols,
#' `+`, `-`, `*`, `/`, `^`, `sqrt`, `exp`, `log`, `softplus`, and `logistic`.
#'
#' `model$param_names`, `model$param_values`, and `model$start_values` retain
#' their legacy role: they contain independent free parameters, in optimization
#' order, on the natural/reporting scale. `model$free_param_names` makes that
#' role explicit. `MCMparameters(model)` returns the complete authoritative
#' table, including fixed and derived parameters, expressions, dependencies,
#' optimizer indices, transformations, and matrix locations.
#'
#' `transform = "positive"` maps an unconstrained optimizer coordinate through
#' softplus and adds `lower` (zero by default). `transform = "bounded"` uses a
#' logistic map between finite `lower` and `upper`. Legacy parameters retain
#' `transform = "identity"` and their existing bound-penalty behavior.
#'
#' Derived expressions are parsed as a restricted abstract syntax tree; they
#' are never evaluated with `eval(parse())`. Unknown symbols, unsupported
#' calls, cycles, invalid starts, and non-finite initial evaluations are errors.
#' Fixed parameters receive SE zero when a covariance is available. Derived
#' SEs and covariances use the delta method from optimizer coordinates.
#'
#' @param model An MCMSEM model.
#' @param name A single parameter name.
#' @param type One of `"free"`, `"fixed"`, or `"derived"`.
#' @param start Natural/reporting-scale starting value for a free parameter.
#' @param value Constant value for a fixed parameter.
#' @param expression One-sided formula for a derived parameter.
#' @param transform Free-parameter transformation: `"identity"`, `"positive"`,
#'   or `"bounded"`.
#' @param lower,upper Natural-scale bounds.
#' @param overwrite Allow replacing an existing definition of the same type.
#' @return A copied model containing the updated parameter graph.
#' @export
#' @examples
#' \dontrun{
#' data <- MCMdatasummary(my_data)
#'
#' # The same signed-gamma cumulant constraint works in the dynamic kernel.
#' dynamic <- MCMmodel(data, kernel = "dynamic")
#' dynamic <- MCMparameter(dynamic, "shape_Earnings", "free", start = 4,
#'                         transform = "positive")
#' dynamic <- MCMparameter(dynamic, "sign_Earnings", "fixed", value = -1)
#' dynamic <- MCMparameter(
#'   dynamic, "tau_Earnings", "derived",
#'   expression = ~ sign_Earnings * 2 / sqrt(shape_Earnings)
#' )
#' dynamic <- MCMparameter(
#'   dynamic, "kappa_Earnings", "derived",
#'   expression = ~ 6 / shape_Earnings
#' )
#'
#' # In a contemporaneous model, k is a raw standardized fourth moment,
#' # hence 3 + excess kurtosis.
#' contemporaneous <- MCMmodel(data, kernel = "contemporaneous")
#' contemporaneous <- MCMparameter(
#'   contemporaneous, "shape_1", "free", start = 4,
#'   transform = "positive"
#' )
#' contemporaneous <- MCMparameter(
#'   contemporaneous, "sk1", "derived", expression = ~ -2 / sqrt(shape_1)
#' )
#' contemporaneous <- MCMparameter(
#'   contemporaneous, "k1", "derived", expression = ~ 3 + 6 / shape_1
#' )
#' }
MCMparameter <- function(model, name, type = c("free", "fixed", "derived"),
                         start, value, expression, transform = NULL,
                         lower = NULL, upper = NULL, overwrite = FALSE) {
  if (!inherits(model, "mcmmodelclass")) {
    stop("`model` must be an MCM model.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name) ||
      !is.na(suppressWarnings(as.numeric(name))) || startsWith(name, "-")) {
    stop("`name` must be one non-numeric parameter label that does not start with `-`.",
         call. = FALSE)
  }
  type <- match.arg(type)
  out <- model$copy()
  .ensure_parameter_graph(out)
  table <- out$parameter_table
  existing <- match(name, table$name)
  is_new <- is.na(existing)
  if (is_new) {
    table <- rbind(table, .parameter_table_row(
      name = name, type = "free", value = 0, start = 0,
      auxiliary = TRUE, matrix_locations = data.frame()
    ))
    existing <- nrow(table)
  } else if (identical(table$type[existing], type) && !isTRUE(overwrite)) {
    stop("Parameter `", name, "` is already defined as `", type,
         "`; use `overwrite = TRUE` to replace that definition.", call. = FALSE)
  }

  current_value <- if (is_new) NA_real_ else table$value[existing]
  table$type[existing] <- type
  table$auxiliary[existing] <- table$auxiliary[existing] ||
    nrow(table$matrix_locations[[existing]]) == 0L
  out$parameter_expressions[[name]] <- NULL

  if (type == "free") {
    start_value <- if (missing(start)) current_value else start
    if (length(start_value) != 1L || !is.finite(start_value)) {
      stop("A free parameter requires one finite `start` value.", call. = FALSE)
    }
    chosen_transform <- if (is.null(transform)) "identity" else match.arg(
      transform, c("identity", "positive", "bounded")
    )
    chosen_lower <- if (is.null(lower)) {
      if (chosen_transform == "positive") 0 else if (is.na(table$lower[existing])) -Inf else table$lower[existing]
    } else lower
    chosen_upper <- if (is.null(upper)) {
      if (is.na(table$upper[existing])) Inf else table$upper[existing]
    } else upper
    if (length(chosen_lower) != 1L || length(chosen_upper) != 1L ||
        is.na(chosen_lower) || is.na(chosen_upper)) {
      stop("`lower` and `upper` must each be one non-missing number.", call. = FALSE)
    }
    table$value[existing] <- table$start[existing] <- as.numeric(start_value)
    table$lower[existing] <- as.numeric(chosen_lower)
    table$upper[existing] <- as.numeric(chosen_upper)
    table$transform[existing] <- chosen_transform
    table$expression[existing] <- ""
    table$dependencies[[existing]] <- character()
  } else if (type == "fixed") {
    fixed_value <- if (missing(value)) current_value else value
    if (length(fixed_value) != 1L || !is.finite(fixed_value)) {
      stop("A fixed parameter requires one finite `value`.", call. = FALSE)
    }
    table$value[existing] <- table$start[existing] <- as.numeric(fixed_value)
    table$lower[existing] <- table$upper[existing] <- NA_real_
    table$transform[existing] <- "none"
    table$expression[existing] <- ""
    table$dependencies[[existing]] <- character()
  } else {
    if (missing(expression)) {
      stop("A derived parameter requires `expression`.", call. = FALSE)
    }
    call <- .parameter_expression_call(expression)
    dependencies <- .parameter_expression_dependencies(call)
    table$transform[existing] <- "none"
    table$lower[existing] <- table$upper[existing] <- NA_real_
    table$expression[existing] <- paste(deparse(call, width.cutoff = 500L),
                                         collapse = " ")
    table$dependencies[[existing]] <- dependencies
    out$parameter_expressions[[name]] <- call
  }
  out$parameter_table <- table
  out$parse()
  out
}

#' Inspect the canonical parameter graph
#'
#' @param object An MCMSEM model or fitted result.
#' @return A data frame with one row per reported parameter.
#' @export
MCMparameters <- function(object) {
  model <- if (inherits(object, "mcmresultclass")) object$model else object
  if (!inherits(model, "mcmmodelclass")) {
    stop("`object` must be an MCM model or result.", call. = FALSE)
  }
  .ensure_parameter_graph(model)
  .parameter_graph_apply_to_model(model)
  .parameter_table_copy(model$parameter_table)
}

.parameter_prune_unused <- function(model) {
  table <- model$parameter_table
  matrix_parameters <- table$name[vapply(table$matrix_locations, function(x) {
    is.data.frame(x) && nrow(x) > 0L
  }, logical(1))]
  needed <- matrix_parameters
  repeat {
    rows <- match(needed, table$name, nomatch = 0L)
    dependencies <- unique(unlist(table$dependencies[rows[rows > 0L]], use.names = FALSE))
    new <- union(needed, dependencies)
    if (setequal(new, needed)) break
    needed <- new
  }
  drop <- table$auxiliary & !(table$name %in% needed)
  if (any(drop)) {
    for (name in table$name[drop]) model$parameter_expressions[[name]] <- NULL
    model$parameter_table <- table[!drop, , drop = FALSE]
  }
  model$parse()
  invisible(model)
}

.parameter_disable_matrix <- function(model, matrix_name) {
  if (!(matrix_name %in% names(model$named_matrices))) return(invisible(model))
  model$named_matrices[[matrix_name]] <- matrix(
    as.character(model$num_matrices[[matrix_name]]),
    nrow = nrow(model$named_matrices[[matrix_name]]),
    ncol = ncol(model$named_matrices[[matrix_name]]),
    dimnames = dimnames(model$named_matrices[[matrix_name]])
  )
  model$parse()
  .parameter_prune_unused(model)
  invisible(model)
}

.parameter_reported_jacobian <- function(model, optimizer_coordinates,
                                         method = "simple") {
  all_names <- model$parameter_table$name
  free_names <- model$param_names
  if (length(optimizer_coordinates) == 0L) {
    return(matrix(0, length(all_names), 0L,
                  dimnames = list(all_names, character())))
  }
  J <- numDeriv::jacobian(
    function(eta) unname(.parameter_values_base(
      model, eta, optimizer_scale = TRUE
    )),
    x = as.numeric(optimizer_coordinates), method = method
  )
  dimnames(J) <- list(all_names, free_names)
  J
}

.parameter_covariance_from_optimizer <- function(model, vcov_optimizer,
                                                 optimizer_coordinates,
                                                 method = "simple") {
  J <- .parameter_reported_jacobian(model, optimizer_coordinates, method)
  V <- J %*% vcov_optimizer %*% t(J)
  V <- (V + t(V)) / 2
  dimnames(V) <- list(model$parameter_table$name, model$parameter_table$name)
  list(
    vcov = V,
    se = stats::setNames(sqrt(pmax(diag(V), 0)), rownames(V)),
    jacobian = J,
    free_vcov = V[model$param_names, model$param_names, drop = FALSE]
  )
}

.parameter_result_table <- function(model, standard_errors = NULL) {
  table <- MCMparameters(model)
  if (is.null(standard_errors)) {
    standard_errors <- stats::setNames(rep(NA_real_, nrow(table)), table$name)
  }
  data.frame(
    parameter = table$name,
    type = table$type,
    estimate = table$value,
    se = unname(standard_errors[table$name]),
    expression = table$expression,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.result_parameter_table <- function(result) {
  stored <- tryCatch(result$parameter_table, error = function(e) data.frame())
  if (is.data.frame(stored) && nrow(stored)) return(stored)
  standard_errors <- if ("se" %in% rownames(result$df)) {
    stats::setNames(as.numeric(result$df["se", ]), colnames(result$df))
  } else NULL
  .parameter_result_table(result$model, standard_errors)
}

# Torch compilation of the same validated graph --------------------------

.parameter_torch_scalar <- function(x, reference) {
  if (inherits(x, "torch_tensor")) return(x)
  torch_tensor(as.numeric(x), device = reference$device, dtype = reference$dtype)
}

.parameter_eval_torch_node <- function(node, values, reference) {
  if (is.numeric(node)) return(.parameter_torch_scalar(node, reference))
  if (is.name(node)) return(values[[as.character(node)]])
  fun <- as.character(node[[1L]])
  args <- lapply(as.list(node)[-1L], .parameter_eval_torch_node,
                 values = values, reference = reference)
  if (fun == "+") return(if (length(args) == 1L) args[[1L]] else args[[1L]] + args[[2L]])
  if (fun == "-") return(if (length(args) == 1L) -args[[1L]] else args[[1L]] - args[[2L]])
  if (fun == "*") return(args[[1L]] * args[[2L]])
  if (fun == "/") return(args[[1L]] / args[[2L]])
  if (fun == "^") return(args[[1L]] ^ args[[2L]])
  if (fun == "sqrt") return(torch_sqrt(args[[1L]]))
  if (fun == "exp") return(torch_exp(args[[1L]]))
  if (fun == "log") return(torch_log(args[[1L]]))
  if (fun == "softplus") {
    x <- args[[1L]]
    return(torch_log1p(torch_exp(-torch_abs(x))) + torch_relu(x))
  }
  if (fun == "logistic") return(torch_sigmoid(args[[1L]]))
  stop("Internal error: unvalidated Torch expression node.", call. = FALSE)
}

.parameter_values_torch <- function(model, optimizer_coordinates) {
  table <- model$parameter_table
  free <- which(table$type == "free")
  if (length(free) != as.numeric(optimizer_coordinates$numel())) {
    stop("Optimizer coordinate tensor has the wrong length.", call. = FALSE)
  }
  values <- list()
  for (j in seq_along(free)) {
    i <- free[j]
    eta <- optimizer_coordinates[j]
    values[[table$name[i]]] <- switch(
      table$transform[i],
      identity = eta,
      positive = table$lower[i] +
        torch_log1p(torch_exp(-torch_abs(eta))) + torch_relu(eta),
      bounded = table$lower[i] + (table$upper[i] - table$lower[i]) *
        torch_sigmoid(eta)
    )
  }
  for (i in which(table$type == "fixed")) {
    values[[table$name[i]]] <- .parameter_torch_scalar(
      table$value[i], optimizer_coordinates
    )
  }
  for (name in .parameter_topological_order(table)) {
    values[[name]] <- .parameter_eval_torch_node(
      model$parameter_expressions[[name]], values, optimizer_coordinates
    )
  }
  values
}

.parameter_torch_matrix <- function(model, values, matrix_name, reference,
                                    named_only = FALSE) {
  numeric_matrix <- model$num_matrices[[matrix_name]]
  base <- if (isTRUE(named_only)) matrix(0, nrow(numeric_matrix), ncol(numeric_matrix)) else numeric_matrix
  table <- model$parameter_table
  for (i in seq_len(nrow(table))) {
    locations <- table$matrix_locations[[i]]
    if (!is.data.frame(locations) || nrow(locations) == 0L) next
    here <- locations$matrix == matrix_name
    if (any(here)) base[locations$index[here]] <- 0
  }
  out <- torch_tensor(base, device = reference$device, dtype = reference$dtype)
  for (i in seq_len(nrow(table))) {
    locations <- table$matrix_locations[[i]]
    if (!is.data.frame(locations) || nrow(locations) == 0L) next
    here <- locations$matrix == matrix_name
    if (!any(here)) next
    map <- matrix(0, nrow(base), ncol(base))
    map[locations$index[here]] <- locations$multiplier[here]
    out <- out + torch_tensor(map, device = reference$device,
                              dtype = reference$dtype) * values[[table$name[i]]]
  }
  out
}

.parameter_torch_reported_vector <- function(model, optimizer_coordinates) {
  values <- .parameter_values_torch(model, optimizer_coordinates)
  torch_stack(values[model$parameter_table$name])
}

# Base-R implied moments compiled from the same graph ---------------------

.contemporaneous_implied_moments_base <- function(model, parameters = NULL,
                                                   optimizer_scale = FALSE,
                                                   use_skewness = TRUE,
                                                   use_kurtosis = TRUE) {
  if (!identical(.model_kernel(model), "contemporaneous")) {
    stop("A contemporaneous MCMSEM model is required.", call. = FALSE)
  }
  out <- model$copy()
  if (!is.null(parameters)) {
    values <- .parameter_values_base(
      out, parameters, optimizer_scale = optimizer_scale
    )
    out$param_values <- unname(values[out$param_names])
  }
  out$inverse_parse()
  A <- out$num_matrices$A
  Fm <- out$num_matrices$Fm
  S <- out$num_matrices$S + 1e-16
  inverse <- solve(diag(nrow(A)) - A)
  M2 <- Fm %*% inverse %*% S %*% t(inverse) %*% t(Fm)
  result <- list(M2 = M2, A = A, Fm = Fm, S = S)
  if (isTRUE(use_skewness)) {
    Sk <- out$num_matrices$Sk
    result$M3 <- Fm %*% inverse %*% Sk %*%
      kronecker(t(inverse), t(inverse)) %*%
      kronecker(t(Fm), t(Fm))
    result$Sk <- Sk
  }
  if (isTRUE(use_kurtosis)) {
    sqrts <- sign(S) * sqrt(abs(S))
    diagonal_S <- S
    diag(diagonal_S) <- 0
    diagonal_S <- all(abs(diagonal_S) < 1e-14)
    K2 <- matrix(1, nrow(out$num_matrices$K), ncol(out$num_matrices$K))
    for (i in seq_len(nrow(S))) {
      coords <- .nd_to_2d_idx(nrow(S), i, i, i, i)
      if (out$named_matrices$K[coords$x, coords$y] != "0") {
        K2[coords$x, coords$y] <- 3
      }
    }
    K_mask <- matrix(1, nrow(out$num_matrices$K), ncol(out$num_matrices$K))
    named_K <- matrix(0, nrow(out$num_matrices$K), ncol(out$num_matrices$K))
    labels <- as.vector(out$named_matrices$K)
    named <- is.na(suppressWarnings(as.numeric(labels)))
    K_mask[named] <- 0
    named_K[named] <- out$num_matrices$K[named]
    K_base <- if (diagonal_S) {
      skron <- kronecker(diag(sqrts), kronecker(diag(sqrts), diag(sqrts)))
      sweep(sqrts %*% out$num_matrices$K1_ref, 2L, skron, `*`)
    } else {
      sqrts %*% out$num_matrices$K1_ref %*%
        kronecker(sqrts, kronecker(sqrts, sqrts))
    }
    K <- K_base * K2 * K_mask + named_K
    result$M4 <- Fm %*% inverse %*% K %*%
      kronecker(kronecker(t(inverse), t(inverse)), t(inverse)) %*%
      kronecker(kronecker(t(Fm), t(Fm)), t(Fm))
    result$K <- K
  }
  result$parameters <- .parameter_values_base(out)
  result
}

#' Calculate model-implied moments in base R
#'
#' Evaluates free, fixed, and derived parameters before constructing either
#' kernel's numerical matrices. This is useful for diagnostics, simulation,
#' and reproducible checks that do not require Torch.
#'
#' @param model An MCMSEM model or fitted result.
#' @param parameters Optional independent free-parameter vector.
#' @param parameter_scale `"reported"` for natural values or `"optimizer"`
#'   for internal coordinates.
#' @param use_skewness,use_kurtosis Include third and fourth moments for a
#'   contemporaneous model. Dynamic models currently require both.
#' @return A named list of implied moment and kernel matrices.
#' @export
MCMimpliedmoments <- function(model, parameters = NULL,
                              parameter_scale = c("reported", "optimizer"),
                              use_skewness = TRUE, use_kurtosis = TRUE) {
  if (inherits(model, "mcmresultclass")) model <- model$model
  if (!inherits(model, "mcmmodelclass")) {
    stop("`model` must be an MCM model or result.", call. = FALSE)
  }
  parameter_scale <- match.arg(parameter_scale)
  if (identical(.model_kernel(model), "dynamic")) {
    if (!isTRUE(use_skewness) || !isTRUE(use_kurtosis)) {
      stop("Dynamic implied moments currently include second through fourth moments.",
           call. = FALSE)
    }
    if (!is.null(parameters) && parameter_scale == "optimizer") {
      values <- .parameter_values_base(model, parameters, optimizer_scale = TRUE)
      parameters <- unname(values[model$param_names])
    }
    result <- .dynamic_implied_moments_base(model, parameters)
    result$parameters <- .parameter_values_base(
      model, if (is.null(parameters)) model$param_values else parameters
    )
    return(result)
  }
  .contemporaneous_implied_moments_base(
    model, parameters, optimizer_scale = parameter_scale == "optimizer",
    use_skewness = use_skewness, use_kurtosis = use_kurtosis
  )
}
