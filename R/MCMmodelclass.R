# Define MCM model class structure
mcmmodelclass <- setRefClass("mcmmodelclass",
                             fields=list(
                               named_matrices="list",
                               num_matrices="list",
                               start_values="data.frame",
                               bounds="data.frame",
                               meta_data="list",
                               param_values="vector",
                               param_names="vector",
                               param_coords="list"
                             ))
# Define MCM model class methods

mcmmodelclass$methods(
    initialize=function(named_matrices=NULL, num_matrices=NULL, start_values=NULL, bounds=NULL, meta_data=NULL,
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
    for (mat in c("A", "Fm", "S")) {
      cat(paste("Matrix", mat, "\n"))
      print(.self$named_matrices[[mat]])
    }

  },
  copy=function() {
    # Create a deepcopy of the model instance
    return(mcmmodelclass(named_matrices=.self$named_matrices, num_matrices=.self$num_matrices, start_values=.self$start_values,
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
        default_u <- .self$meta_data$bound_default[["U"]][[sub("^([[:alpha:]]*).*", "\\1", gsub("l", "",param))]]
        default_l <- .self$meta_data$bound_default[["L"]][[sub("^([[:alpha:]]*).*", "\\1", gsub("l", "",param))]]
        new_col <- data.frame(rbind(default_l, default_u))
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
            .self$bounds[row, col] <- .self$meta_data$bound_default[[row]][[sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", col))]]
          }
        }
      }
    }
    .self$start_values <- as.data.frame(t(.self$param_values))
    colnames(.self$start_values) <- .self$param_names
    rownames(.self$start_values) <- "start"
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

print.mcmmodelclass <- function(model, matrix=NULL) {
  if (is.null(matrix)) {
    model$show()
  } else {
    if (matrix %in% c("A", "Fm", "S", "Sk", "K")) {
      print(model$named_matrices[[matrix]])
    } else {
      stop('matrix argument should be one of ("A", "Fm", "S", "Sk", "K")')
    }
  }
}