.string_parse <- function(x) {
  # Parse a string of "1*1+0*2+sqrt(a)*sqrt(a)+a*sqrt(0)+1*a*b+3*2" to "1+a+a*b+6"
  out <- c()
  for (i in strsplit(x, "\\+")[[1]]) {
    if ("0" %in% strsplit(i, "\\*")[[1]]) {
      next
    } else if ("1" %in% strsplit(i, "\\*")[[1]]) {
      i <- paste0(strsplit(i, "\\*")[[1]][strsplit(i, "\\*")[[1]] != "1"], collapse="*")
    }
    if (grepl("sqrt\\(", i)) {
      # Handle square roots
      pars <- gsub("\\)", "", gsub("sqrt\\(", "", regmatches(i, gregexpr("sqrt\\(.*?\\)", i))[[1]]))
      nosqrtstr <- paste0(strsplit(i, "\\*")[[1]][!startsWith(strsplit(i, "\\*")[[1]], "sqrt(")], collapse="*")
      if ("0" %in% pars) {
        next
      }
      if (length(pars) > 2 ) {
        j <- c()
        for (par in unique(pars)) {
          if (table(pars)[par] == 2) {
            j <- c(j, par)
          } else {
            if (par != "1") {
              j <- c(j, paste0("sqrt(",par,")"))
            }
          }
        }
        if (nosqrtstr == "") {
          i <- paste0(j, collapse="*")
        } else {
          i <- paste0(nosqrtstr, "*", paste0(j, collapse="*"))
        }
      } else if (length(pars) == 2) {
        if (pars[1] == pars[2]) {
          sqrtout <- pars[1]
        } else if (pars[1] == "1") {
          if (pars[2] == 1) {
            sqrtout <- "1"
          } else {
            sqrtout <- paste0("sqrt(", pars[2], ")")
          }
        } else if (pars[2] == "1") {
          sqrtout <- paste0("sqrt(", pars[1], ")")
        } else {
          out <- c(out, i)
          next
        }
        if (nosqrtstr == "") {
          i <- sqrtout
        } else {
          i <- paste0(nosqrtstr, "*", sqrtout)
        }
      } else {
        # If only one side is square root:
        # Replace sqrt(1) with 1 and sqrt(0) with 0, and continue as normal
        if (pars == "0") {
          i <- gsub("sqrt\\(0\\)", "0", i)
        } else if (pars == "1") {
          i <- gsub("sqrt\\(1\\)", "1", i)
        }
      }
    }
    if (all(is.finite(suppressWarnings(as.numeric(strsplit(i, "\\*")[[1]]))))) {
      # if i consists of only numbers:
      i_num <- as.numeric(strsplit(i, "\\*")[[1]])
      if (all(i_num != 0)) {
        res <- 1
        for (j in seq_along(i_num)) {
          res <- res * i_num[j]
        }
        out <- c(out, as.character(res))
      }
    } else if (!(startsWith(i, "0*") | endsWith(i, "*0"))) {
      if (startsWith(i, "1*")) {
        out <- c(out, gsub("1\\*", "", i))
      } else if (endsWith(i, "*1")) {
        out <- c(out, gsub("\\*1", "", i))
      } else {
        out <- c(out, i)
      }
    }
  }
  if (length(out) > 0) {
    return(paste0(out, collapse="+"))
  } else {
    return("0")
  }
}

.string_sqrt <- function(a) {
  matrix(paste0("sqrt(", a, ")"), nrow=nrow(a), ncol=ncol(a))
}

.string_mulvec <- function(a, b) {
  # string-multiply 2 vectors, i.e. c(1, a) c(b, 2) to c(1*b, a*2)
  out <- c()
  for (i in seq_along(a)) {
    if (grepl("\\+", a[i])) {
      if (grepl("\\+", b[i])) {
        outi <- c()
        for (j in strsplit(a[i], "\\+")[[1]]) {
          for (k in strsplit(b[i], "\\+")[[1]]) {
            outi <- c(outi, paste0(j, "*",k))
          }
        }
        out <- c(out, paste0(outi, collapse="+"))
      } else {
        outi <- c()
        for (j in strsplit(a[i], "\\+")[[1]]) {
          outi <- c(outi, paste0(j, "*", b[i]))
        }
        out <- c(out, paste0(outi, collapse="+"))
      }
    } else if (grepl("\\+", b[i])) {
      outi <- c()
      for (j in strsplit(b[i], "\\+")[[1]]) {
        outi <- c(outi, paste0(a[i], "*", j))
      }
      out <- c(out, paste0(outi, collapse="+"))
    } else {
      out <- c(out, paste0(a[i], "*", b[i]))
    }
  }
  return(out)
}

.string_matmul <- function(a, b) {
  # Matrix mulplication %*% with strings, using the two previously defined functions
  out <- matrix("0", nrow=nrow(a), ncol=ncol(b))
  for (i in seq_len(nrow(a))) {
    for (j in seq_len(ncol(b))) {
      out[i, j] <- .string_parse(paste0(.string_mulvec(a[i, ], b[, j]), collapse="+"))
    }
  }
  return(out)
}

.string_kron <- function(a, b) {
  # Kronecker product %x% of string matrices using previously defined funcitons
  out <- matrix("0", nrow=nrow(a)*nrow(b), ncol=ncol(a)*ncol(b))
  for (i in seq_len(nrow(b))) {
    for (j in seq_len(ncol(b))) {
      curmatrix <- matrix(.string_mulvec(rep(b[i, j], nrow(a)*ncol(a)), as.vector(a)), nrow=nrow(a), ncol=ncol(a))
      for (ii in seq_len(nrow(curmatrix))) {
        for (jj in seq_len(ncol(curmatrix))) {
          curmatrix[ii, jj] <- .string_parse(curmatrix[ii, jj])
        }
      }
      out[((i-1)*nrow(a)+1):(i*nrow(a)), ((j-1)*ncol(a)+1):(j*ncol(a))] <- curmatrix
    }
  }
  return(out)
}

.string_formshorten <- function(x) {
  # Shorten such that a*b+a*b becomes 2*(a*b)
  tab <- table(strsplit(x, "\\+")[[1]])
  out <- c()
  for (i in names(table(strsplit(x, "\\+")[[1]]))) {
    if (any(table(strsplit(i, "\\*")) > 1)) {
      tab2 <- table(strsplit(i, "\\*"))
      lab <- c()
      for (j in names(tab2)) {
        if (unname(tab2[j]) > 1) {
          lab <- c(lab, paste0(j,"^",unname(tab2[j])))
        } else {
          lab <- c(lab, j)
        }
      }
      lab <- paste0(lab, collapse="*")
    } else {
      lab <- i
    }
    if (unname(tab[i]) > 1) {
      if (grepl("\\+", lab)) {
        out <- c(out, paste0(unname(tab[i]),"*(",lab,")"))
      } else {
        out <- c(out, paste0(unname(tab[i]),"*",lab))
      }
    } else {
      out <- c(out, lab)
    }
  }
  return(paste0(out,collapse="+"))
}
