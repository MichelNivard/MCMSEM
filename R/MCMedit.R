MCMedit <- function(model, pointer, name, value) {
  x <- model$copy()
  if (pointer %in% names(x$num_matrices)) {
    # For modifying parameters:
    # Verify integrity of named parameters, they should be a, b, s, sk or k if all but letters are dropped
    if (is.character(value)) {
      if (length(value) > 1) {
        for (i in value) {
          i_pos <- gsub("-", "", i)
          if (!(sub("^([[:alpha:]]*).*", "\\1", i_pos) %in% c("a", "s", "b", "sk", "k", "fm")))
            stop("Named parameters must only contain one of [a, s, b, sk, k, fm] and numbers or other symbols")
        }
      } else {
        value_pos <- gsub("-", "", value)
        if (!(sub("^([[:alpha:]]*).*", "\\1", value_pos) %in% c("a", "s", "b", "sk", "k", "fm")))
          stop("Named parameters must only contain one of [a, b, s, sk, k, fm] and numbers or other symbols")
      }
    }
    if (length(name) > 1) {
      if (is.character(value)) {
        x$named_matrices[[pointer]][name[1], name[2]] <- value
      } else {
        x$num_matrices[[pointer]][name[1], name[2]] <- value
        old_name <- x$named_matrices[[pointer]][name[1], name[2]]
        x$named_matrices[[pointer]][name[1], name[2]] <- as.character(value)
        x$bounds[, old_name] <- NULL # Since paramter is set to a constant: Remove bounds
      }
    } else {
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