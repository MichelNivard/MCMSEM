MCMedit <- function(model, pointer, name, value) {
  x <- model$copy()
  if (pointer %in% names(x$num_matrices)) {
    # For modifying parameters:
    # Verify integrity of named parameters, they should be a, b, s, sk or k if all but letters are dropped
    if (is.character(value) & all(is.na(suppressWarnings(as.numeric(value))))) {
      if (length(value) > 1) {
        for (i in value) {
          i_pos <- gsub("-", "", i)
          if (!(sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", i_pos)) %in% c("a", "s", "b", "sk", "k", "fm")))
            stop("Named parameters must only contain one of [a, s, b, sk, k, fm] and numbers or other symbols")
        }
      } else {
        value_pos <- gsub("-", "", value)
        if (!(sub("^([[:alpha:]]*).*", "\\1", gsub("l", "", value_pos)) %in% c("a", "s", "b", "sk", "k", "fm")))
          stop("Named parameters must only contain one of [a, b, s, sk, k, fm] and numbers or other symbols")
      }
    }
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
          }
        } else {
          for (i in seq_along(name[[1]])) {
            x$num_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- value[i]
            old_name <- x$named_matrices[[pointer]][name[[1]][i], name[[2]][i]]
            x$named_matrices[[pointer]][name[[1]][i], name[[2]][i]] <- as.character(value[i])
            x$bounds[, old_name] <- NULL # Since paramter is set to a constant: Remove bounds
          }
        }
      } else {
        if (is.character(value)) {
          x$named_matrices[[pointer]][name[1], name[2]] <- value
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
            x$bounds[, old_name] <- NULL # Since paramter is set to a constant: Remove bounds
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
      }
      if (is.character(value)) {
        x$named_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- value
      } else if (is.numeric(value)) {
        x$num_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- value
        x$named_matrices[[pointer]][x$named_matrices[[pointer]] == name] <- as.character(value)
      }
    }
  } else if (pointer %in% c("bound", "ubound", "lbound")){
    # For modifying bounds
    if (name %in% c("a", "b", "s", "sk", "k", "fm")) {
      col_to_change <- which(sub("^([[:alpha:]]*).*", "\\1", colnames(x$bounds)) == name)
    } else if (name %in% colnames(x$bounds)) {
      col_to_change <- which(colnames(x$bounds) == name)
    }
    row_to_change <- list(bound=c(1, 2), lbound=1, ubound=2)[[pointer]]
    x$bounds[row_to_change, col_to_change] <- value
  } else if (pointer == "start") {
    if ((!(name %in% x$param_names)) & !(name %in% c("a", "b", "s", "sk", "k", "fm")) ) {
      stop(paste0("Parameter ", name, " not found"))
    } else if (name %in% c("a", "b", "s", "sk", "k", "fm")) {
      cols_to_change <- which(sub("^([[:alpha:]]*).*", "\\1", colnames(x$start_values)) == name)
      x$start_values["start", cols_to_change] <- value
      x$param_values[cols_to_change] <- value
    } else {
      x$start_values[name] <- value
      x$param_values[x$param_names == name] <- value
    }
    x$inverse_parse()
  } else {
    stop("Second input argument not recognized")
  }
  x$parse()
  return(x)
}