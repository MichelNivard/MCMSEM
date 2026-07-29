MCMedit <- function(model, pointer, name, value) {
  default_starts <- list(a=0.2, b=0, s=0, sk=0, k=0)
  matrix_start <- function(label, current) {
    prefix <- sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", gsub("-", "", label)))
    proposed <- default_starts[[prefix]]
    if (is.null(proposed) || length(proposed) == 0L || !is.finite(proposed)) {
      proposed <- if (length(current) && is.finite(current)) current else 0
    }
    proposed
  }
  x <- model$copy()
  dynamic_model <- identical(.model_kernel(x), "dynamic")
  if (pointer %in% names(x$num_matrices)) {
    if (length(name) > 1) {
      if (is.list(name)) {
        if (length(name[[1]]) == 1) {
          name[[1]] <- rep(name[[1]], length(name[[2]]))
        } else if (length(name[[2]]) == 1) {
          name[[2]] <- rep(name[[2]], length(name[[1]]))
        }
        if (length(value) == 1) {
          value <- rep(value, length(name[[1]]))
        } else if (length(value) != length(name[[1]])) {
          stop("Length of provided names/coordinates is not equal to length of replacements.")
        }

        if (is.character(value)) {
          for (i in seq_along(name[[1]])) {
            x$named_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- value[i]
            if (!dynamic_model) {
              current <- x$num_matrices[[pointer]][name[[1]][i], name[[2]][i]]
              x$num_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- matrix_start(value[i], current)
            }
          }
        } else {
          for (i in seq_along(name[[1]])) {
            x$num_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- value[i]
            old_name <- x$named_matrices[[pointer]][name[[1]][i], name[[2]][i]]
            x$named_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- as.character(value[i])
            if (old_name %in% colnames(x$bounds)) x$bounds[, old_name] <- NULL
          }
        }
      } else {
        if (is.character(value)) {
          x$named_matrices[[pointer]][name[1], name[2]] <- value
          if (!dynamic_model) {
            current <- x$num_matrices[[pointer]][name[1], name[2]]
            x$num_matrices[[pointer]][name[1], name[2]] <- matrix_start(value, current)
          }
        } else {
          if (is.character(name)) {
            # MCMedit(model, "A", c("b1_2", "b1_3", "b1_4"), 0)
            for (i in name) {
              if (is.character(value)) {
                x$named_matrices[[pointer]][gsub("-", "", x$named_matrices[[pointer]]) == i] <- value
              } else if (is.numeric(value)) {
                x$num_matrices[[pointer]][gsub("-", "", x$named_matrices[[pointer]])== i] <- value
                x$named_matrices[[pointer]][gsub("-", "", x$named_matrices[[pointer]]) == i] <- as.character(value)
              }
            }
          } else {
            # MCMedit(model, "A", c(1, 2), 0)
            x$num_matrices[[pointer]][name[1], name[2]] <- value
            old_name <- x$named_matrices[[pointer]][name[1], name[2]]
            x$named_matrices[[pointer]][name[1], name[2]] <- as.character(value)
            if (old_name %in% colnames(x$bounds)) x$bounds[, old_name] <- NULL
          }
        }
      }
    } else {
      if (name %in% c("a", "s", "b", "sk", "k", "fm")) {
        idx <- which(startsWith(gsub("-", "", x$named_matrices[[pointer]]), name))
        if (is.character(value)) {
          x$named_matrices[[pointer]][idx] <- value
        } else if (is.numeric(value)) {
          x$num_matrices[[pointer]][idx] <- value
          x$named_matrices[[pointer]][idx] <- as.character(value)
        }
      } else {
        if (is.character(value)) {
          x$named_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- value
        } else if (is.numeric(value)) {
          x$num_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- value
          x$named_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- as.character(value)
        }
      }
    }
  } else if (pointer %in% c("bound", "ubound", "lbound")){
    # For modifying bounds
    if (dynamic_model && all(name %in% c("B", "Tau", "Kappa", "L_G"))) {
      group_by_parameter <- vapply(x$param_coords, `[[`, character(1), 1)
      col_to_change <- which(group_by_parameter %in% name)
    } else if (all(name %in% c("a", "b", "s", "sk", "k", "fm"))) {
      colsub <- sub("^([[:alpha:]]*).*", "\\1", colnames(x$bounds))
      colsub[startsWith(colnames(x$bounds), "sk")] <- "sk"
      col_to_change <- which(colsub %in% name)
    } else if (all(name %in% colnames(x$bounds))) {
      col_to_change <- which(colnames(x$bounds) == name)
    } else {
      defined <- intersect(name, x$parameter_table$name)
      if (length(defined)) {
        stop("Bounds can only be edited for free parameters; `", defined[1],
             "` is ", x$parameter_table$type[match(defined[1], x$parameter_table$name)],
             ".", call. = FALSE)
      }
      stop("Parameter not found in model bounds: ", paste(name, collapse = ", "),
           call. = FALSE)
    }
    row_to_change <- list(bound=c(1, 2), lbound=1, ubound=2)[[pointer]]
    x$bounds[row_to_change, col_to_change] <- value
  } else if (pointer == "start") {
    if (name == "all") {
      if (!(length(value) %in% c(x$start_values$getncol(), 1))) {
        stop(paste0("Value should either be of length 1, or ", x$start_values$getncol()))
      } else {
        x$start_values$set_all(value)
        x$param_values <- value
      }
    } else if (dynamic_model && all(name %in% c("B", "Tau", "Kappa", "L_G"))) {
      group_by_parameter <- vapply(x$param_coords, `[[`, character(1), 1)
      cols_to_change <- which(group_by_parameter %in% name)
      x$start_values$set("start", cols_to_change, value)
      x$param_values[cols_to_change] <- value
    } else if (all((!(name %in% x$param_names)) & !(name %in% c("a", "b", "s", "sk", "k", "fm")) )) {
      graph_match <- match(name, x$parameter_table$name)
      if (!is.na(graph_match)) {
        stop("Starting values can only be edited for free parameters; `", name,
             "` is ", x$parameter_table$type[graph_match], ".", call. = FALSE)
      }
      stop(paste0("Parameter ", name, " not found"))
    } else if (all(name %in% c("a", "b", "s", "sk", "k", "fm"))) {
      cols_to_change <- which(sub("^([[:alpha:]]*).*", "\\1", x$start_values$getcolnames()) == name)
      x$start_values$set("start", cols_to_change, value)
      x$param_values[cols_to_change] <- value
    } else {
      x$start_values$set("start", x$param_names == name, value)
      x$param_values[x$param_names == name] <- value
    }
    x$inverse_parse()
  } else {
    stop("Second input argument not recognized")
  }
  x$parse()
  return(x)
}
