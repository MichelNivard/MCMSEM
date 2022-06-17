# Wrapper function to create mcmmodel instance
MCMmodel <- function(data, n_confounding=1, constrained_a=TRUE, scale_data=TRUE) {
  # TODO: Expand checks on how many confoundings can/should be used with or without constrained a depending on input data
  ## Make matrices with names
  if (nrow(data) < 1000)
    stop("Currently only a dataframe with at least 1000 rows is supported.")
  if (ncol(data) < 2) {
    stop("At least two columns in data are required.")
  } else if (ncol(data) == 2) {
    if (n_confounding > 1) {
      stop("Only one confounding can be used with two variables.")
    } else if (!(constrained_a)) {
      stop("Constrained a has to be used with two variables.")
    }
  }
  if (n_confounding > ncol(data)) {
    stop("Cannot use more confounding factors than phenotypes present in the data")
  }
  # Scale data
  data_scaled <- apply(data, 2, scale)
  data_was_scaled <- all(round(data, 2) == round(data_scaled, 2))
  if ((scale_data) & !(data_was_scaled)) {
    data <- data_scaled
  }

  data <- as.matrix(data)

  n_f <- n_confounding
  n_p <- ncol(data)
  if (constrained_a) {
    a_names <- c()
    for (i in 1:n_f) {
      a_names <- c(a_names, rep(paste0("a", i), n_p))
    }
  } else {
    a_names <- c()
    for (i in 1:n_f) {
      for (j in 1:n_p) {
        a_names <- c(a_names, paste0("a", i, "_", j))
      }
    }
  }
  b_names <- c()
  for (i in 1:n_p) {
    js <- 1:n_p
    for (j in js[js != i]) {
      b_names <- c(b_names, paste0("b",j,"_", i))
    }
  }

  par_names <- list(
    a=a_names,
    b=b_names,
    s=paste0("s", 1:n_p),
    sk=paste0("sk", 1:n_p),
    k=paste0("k", 1:n_p)
  )

  named_matrices <- .gen_matrices(par_names, n_p, n_f, base_value="0")

  ## Make matrices with starting values
  M3.obs <- M3.MM(data)
  M4.obs <- M4.MM(data)
  sk_starts <- c()
  k_starts <- c()
  for (i in 1:n_p) {
    sk_starts <- c(sk_starts, M3.obs[i, i + (i-1)*n_p])
    k_starts <- c(k_starts, M4.obs[i, i + (i-1)*n_p + (i-1)*n_p^2])
  }
  par <- list(
    a=rep(0.2, length(par_names[['a']])),
    b=rep(0.0, length(par_names[['b']])),
    s=rep(1, length(par_names[['s']])),
    sk=sk_starts,
    k=k_starts
  )
  # Lower and upper bounds, defaults are provided here, values are filled in mcmmodelclass$parse()
  L <- rep(NA, length(unname(unlist(par_names))))
  U <- rep(NA, length(unname(unlist(par_names))))
  bounds <- as.data.frame(rbind(L, U))
  colnames(bounds) <- unname(unlist(par_names))
  bound_defaults <- list(L=list(a=-1, b=-.5, s=0.01, sk=-5, k=0, fm=-1),
                         U=list(a=1, b=1, s=2, sk=18, k=100, fm=1))
  num_matrices <- .gen_matrices(par, n_p, n_f, base_value=0)
  start_values <- as.data.frame(t(unlist(par)))
  colnames(start_values) <- unlist(par_names)
  rownames(start_values) <- "start"
  if (n_confounding > 1) {
    for (n_fi in 1:(n_confounding-1)) {
      num_matrices[["A"]][(n_confounding+1):(n_confounding+n_fi), n_fi+1] <- 0.0
      named_matrices[["A"]][(n_confounding+1):(n_confounding+n_fi), n_fi+1] <- "0"
    }
  }
  model <- mcmmodelclass(named_matrices=named_matrices,
                       num_matrices=num_matrices,
                       start_values=start_values,
                       bounds=bounds,
                       meta_data=list(n_phenotypes=n_p, n_confounding=n_f, bound_defaults=bound_defaults,
                                      data_was_scaled=data_was_scaled, scale_data=scale_data))

  return(model)
}