MCMparseK <- function(model, sort=TRUE, shorten=TRUE, add_one=TRUE, print=FALSE) {
  S <- model$named_matrices$S
  if (ncol(S) > 11) {
    cat("This is probably going to take a while... \nFeel free to cancel this operation, or grab yourself a drink and a snack (or a few).\n")
  }
  # v1-style: K[ round(M4.MM(cbind(rnorm(2000), rnorm(2000), rnorm(2000)))) != 0] <- 1
  K <- matrix(as.character(model$num_matrices$K1_ref+1-1), nrow=nrow(model$num_matrices$K1_ref), ncol=ncol(model$num_matrices$K1_ref))
  # v1-style: K <- sqrt(S) %*% K %*% (sqrt(S) %x% sqrt(S) %x% sqrt(S))
  Kstr <- .string_matmul(.string_matmul(.string_sqrt(S), K), .string_kron(.string_kron(.string_sqrt(S), .string_sqrt(S)), .string_sqrt(S)))
  if (sort | shorten) {
    for (idx in which(Kstr != "0")) {
      if (sort) {
        # Sort parameters such that s2*s1 becomes s1*s2
        Kstr[idx] <- do.call(paste0, list(lapply(lapply(lapply(strsplit(strsplit(Kstr[idx], "\\+")[[1]], "\\*"), sort), paste0, collapse="*"), sort), collapse="+"))
      }
      if (shorten) {
        # Shorten such that a*b+a*b becomes 2*(a*b)
        Kstr[idx] <- .string_formshorten(Kstr[idx])
      }
    }
  }
  if (add_one) {
    # Add 1* to single parameters
    Kstr[!(Kstr %in% c("0", "1") | grepl("\\*", Kstr)  | grepl("\\+", Kstr))] <- paste0("1*", Kstr[!(Kstr %in% c("0", "1") | grepl("\\*", Kstr) | grepl("\\+", Kstr))])
  }
  # Add actual parameters from the K matrix
  Kstr[model$named_matrices$K != "0"] <- model$named_matrices$K[model$named_matrices$K != "0"]
  if (print) {
    # Print parsed K formatted such that it is as readable as it gets
    formattedK <- apply(Kstr, c(1,2), function(x, maxnchar) {if (nchar(x) < maxnchar) {return(paste0(x, paste0(rep(' ', maxnchar-nchar(x)), collapse='')))} else {return(x)}}, maxnchar=max(nchar(Kstr)))
    print(formattedK, quote=FALSE)
    return(invisible(NULL))
  } else {
    return(Kstr)
  }
}
