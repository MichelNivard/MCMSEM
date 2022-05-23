mcmedit <- function(x, y, z, value) {
  if (y %in% names(x$num_matrices)) {
    # For modifying parameters:
    # Verify integrity of named parameters, they should be a, b, s, sk or k if all but letters are dropped
    if (is.character(value)) {
      if (length(value) > 1) {
        for (i in value) {
          if (!(sub("^([[:alpha:]]*).*", "\\1", i) %in% c("a", "b", "sk", "k", "fm")))
            stop("Named parameters must only contain one of [a, b, sk, k, fm] and numbers or other symbols")
        }
      } else {
        if (!(sub("^([[:alpha:]]*).*", "\\1", value) %in% c("a", "b", "sk", "k", "fm")))
          stop("Named parameters must only contain one of [a, b, sk, k, fm] and numbers or other symbols")
      }
    }
    if (length(z) > 1) {
      if (is.character(value)) {
        x$named_matrices[[y]][z[1], z[2]] <- value
      } else {
        x$num_matrices[[y]][z[1], z[2]] <- value
        old_name <- x$named_matrices[[y]][z[1], z[2]]
        x$named_matrices[[y]][z[1], z[2]] <- as.character(value)
        x$bounds[, old_name] <- NULL # Since paramter is set to a constant: Remove bounds
      }
    } else {
      if (is.character(value)) {
        x$named_matrices[[y]][x$named_matrices[[y]] == z] <- value
      } else if (is.numeric(value)) {
        x$num_matrices[[y]][x$named_matrices[[y]] == z] <- value
        x$named_matrices[[y]][x$named_matrices[[y]] == z] <- as.character(value)
      }
    }
  } else if (y %in% c("bound", "ubound", "lbound")){
    # For modifying bounds
    if (z %in% c("a", "b", "s", "sk", "k", "fm")) {
      col_to_change <- which(sub("^([[:alpha:]]*).*", "\\1", colnames(x$bounds)) == z)
    } else if (z %in% colnames(x$bounds)) {
      col_to_change <- which(colnames(x$bounds) == z)
    }
    row_to_change <- list(bound=c(1, 2), lbound=1, ubound=2)[[y]]
    x$bounds[row_to_change, col_to_change] <- value
  } else if (y == "start") {
    if (!(z %in% x$param_names))
      stop(paste0("Parameter ", z, " not found"))
    x$param_values[x$param_names == z] <- value
  } else {
    stop("Second input argument not recognized")
  }
  x$parse()
  return(invisible(NULL))
}