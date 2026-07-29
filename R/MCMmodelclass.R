# Custom start values class as a workaround to prevent editing through model$start_values["start", "a"] <- 1.0
mcmstartvaluesclass <- setRefClass("mcmstartvaluesclass",
                                   fields=list(.df="data.frame"))
mcmstartvaluesclass$methods(
  initialize=function(data=data.frame(NULL)) {.self$.df <- data},
  copy=function() {return(mcmstartvaluesclass(.self$.df))},
  show=function() {print(.self$.df)},
  get=function(x, y) {return(.self$.df[x, y])},
  getncol=function() {return(ncol(.self$.df))},
  set=function(x, y, z) {.self$.df[x, y] <- z},
  set_all=function(x) {.self$.df["start", ] <- x},
  getcolnames=function() {return(colnames(.self$.df))},
  drop=function(idx) {.self$.df <- .self$.df[, -idx]}
)

.model_default_bounds <- function(model, parameter) {
  idx <- match(parameter, model$param_names)
  coord <- model$param_coords[[idx]]
  matrix_name <- coord[[1]]
  if (identical(.model_kernel(model), "dynamic")) {
    if (matrix_name == "B") {
      rc <- lapply(coord[[2]], .r_1to2d_idx,
                   nrows = nrow(model$num_matrices$B))
      diagonal <- all(vapply(rc, function(x) x[1] == x[2], logical(1)))
      return(c(L = if (diagonal) 0 else -0.98, U = 0.98))
    }
    if (matrix_name == "Tau") return(c(L = -100, U = 100))
    if (matrix_name == "Kappa") return(c(L = -1.99, U = 100))
    if (matrix_name == "L_G") {
      rc <- lapply(coord[[2]], .r_1to2d_idx,
                   nrows = nrow(model$num_matrices$L_G))
      diagonal <- all(vapply(rc, function(x) x[1] == x[2], logical(1)))
      return(c(
        L = if (diagonal) log(1e-6) else -100,
        U = if (diagonal) log(1e2) else 100
      ))
    }
  }
  defaults <- model$meta_data$bound_default
  if (is.null(defaults)) defaults <- model$meta_data$bound_defaults
  prefix <- sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", parameter))
  c(L = defaults$L[[prefix]], U = defaults$U[[prefix]])
}

# Define MCM model class structure
mcmmodelclass <- setRefClass("mcmmodelclass",
                             fields=list(
                               named_matrices="list",
                               num_matrices="list",
                               start_values="mcmstartvaluesclass",
                               bounds="data.frame",
                               meta_data="list",
                               param_values="vector",
                               param_names="vector",
                               param_coords="list"
                             ))

# Define MCM model class methods
mcmmodelclass$methods(
    initialize=function(named_matrices=NULL, num_matrices=NULL, start_values=mcmstartvaluesclass(), bounds=NULL, meta_data=NULL,
                        param_values=c(0), param_names=c(""), param_coords=list()){
      if (!(is.null(named_matrices))) {
        # This is executed upon initialization, required to force parse upon initialization of class instance
        .self$named_matrices <- named_matrices
        .self$num_matrices <- num_matrices
        .self$start_values <- start_values
        .self$bounds <- bounds
        .self$meta_data <- meta_data
        .self$param_values <- param_values
        .self$param_names <- param_names
        .self$param_coords <- param_coords
        if (all(param_names == c(""))) {
          # When a new class is made: parse, if a copy is made parsing is not necessary
          .self$parse()
        }
      }
  },
  show=function() {
    # This is just what shows when you run 'mcmmodelinstance' in command prompt
    kernel <- .model_kernel(.self)
    cat("Kernel: ", .kernel_label(kernel), "\n", sep = "")
    matrices <- if (identical(kernel, "dynamic")) {
      c("B", "Tau", "Kappa", "L_G", "D2")
    } else {
      c("A", "Fm", "S")
    }
    for (mat in matrices) {
      cat(paste("Matrix", mat, "\n"))
      print(.self$named_matrices[[mat]])
    }

  },
  copy=function() {
    # Create a deepcopy of the model instance
    return(mcmmodelclass(named_matrices=.self$named_matrices, num_matrices=.self$num_matrices, start_values=.self$start_values$copy(),
                         bounds=.self$bounds, meta_data=.self$meta_data, param_values=.self$param_values, param_names=.self$param_names,
                         param_coords=.self$param_coords))
  },
  parse=function() {
    # Parse named matrices to param_values, param_coords and bounds (to make it easy to access these during optimization)
    # Update param_values and param_coords
    .self$param_values <- c(0)  # vector cannot be empty due to it being set as vector in class fields
    .self$param_names <- c("")
    .self$param_coords <- list()
    for (mat in names(.self$named_matrices)) {
      # Ensure paramters starting with - are not recognized as unique parameter and store in vector for easier indexing
      neg_indices <- startsWith(as.vector(.self$named_matrices[[mat]]), "-") # Which parameter are labeled -
      non_neg_vec <- gsub("-", "", as.vector(.self$named_matrices[[mat]]))  # Vectorized named matrix without "-"labels
      unique_params <- unique(non_neg_vec)
      unique_params <- unique_params[is.na(suppressWarnings(as.numeric(unique_params)))]
      for (param in unique_params) {
        current_coords <- which(non_neg_vec == param)
        # neg_indices[current_coords]*-2+1 produces multipliers for parameters labeled -, e.g. "-a1":
        #  if the matrix contains "-a1, a1" in whatever place, neg_indices will be TRUE FALSE there respectively
        #  c(1, 0) * -2 + 1 produces (-1, 1), i.e. multipliers for these paramter values.
        .self$param_coords <- append(.self$param_coords, list(list(mat, current_coords, neg_indices[current_coords]*-2+1)))
        .self$param_names <- c(.self$param_names, param)
        .self$param_values <- c(.self$param_values, .self$num_matrices[[mat]][current_coords][1])
      }
    }
    .self$param_names <- .self$param_names[2:length(.self$param_names)]
    .self$param_values <- .self$param_values[2:length(.self$param_values)]
    # Update bounds
    for (param in .self$param_names) {
      if (!(param %in% colnames(.self$bounds))) {
        default <- .model_default_bounds(.self, param)
        new_col <- data.frame(rbind(default[["L"]], default[["U"]]))
        colnames(new_col) <- param
        .self$bounds <- cbind(.self$bounds, new_col)
      }
    }
    for (col in colnames(.self$bounds)) {
      if (!(col %in% .self$param_names)) {
        .self$bounds[, col] <- NULL
      } else {
        if (any(is.na(.self$bounds[, col]))) {
          for (row in rownames(.self$bounds)[which(is.na(.self$bounds[, col]))]) {
            .self$bounds[row, col] <- .model_default_bounds(.self, col)[[row]]
          }
        }
      }
    }
    starts <- as.data.frame(t(.self$param_values))
    colnames(starts) <- .self$param_names
    rownames(starts) <- "start"
    .self$start_values <- mcmstartvaluesclass(starts)
    .self$bounds <- .self$bounds[.self$param_names]
  },
  inverse_parse=function() {
    for (n_par in seq_along(.self$param_names)) {
      par_name <- .self$param_names[n_par]
      for (mat in names(.self$named_matrices)) {
        .self$num_matrices[[mat]][.self$named_matrices[[mat]] == par_name] <- .self$param_values[n_par]
      }
    }
  }
)

print.mcmmodelclass <- function(x, matrix=NULL, ...) {
  model <- x
  if (is.null(matrix)) {
    model$show()
  } else {
    available <- names(model$named_matrices)
    if (matrix %in% setdiff(available, c("K"))) {
      print(model$named_matrices[[matrix]])
    } else if (matrix == "K" && "K" %in% available) {
      print(model$named_matrices[[matrix]])
      warning(paste0("This matrix only contains free K parameters, the K matrix used for optimization is additionally a product of the S matrix. To see this full product, run MCMparseK(",deparse(substitute(model)),")"))
    } else {
      stop(
        "`matrix` must be one of: ",
        paste(available, collapse = ", "),
        call. = FALSE
      )
    }
  }
}

summary.mcmmodelclass <- function(object, ...) {
  object$show()
  dof <- MCMdegreesoffreedom(object)
  cat("Free parameters: ", dof$n_parameters, "\n", sep = "")
  cat("Unique moments: ", dof$n_moments, "\n", sep = "")
  cat("Nominal degrees of freedom: ", dof$df, "\n", sep = "")
  invisible(object)
}
