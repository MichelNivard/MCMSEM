# Wrapper function to create mcmmodel instance
MCMmodel <- function(data, n_latent=1, constrained_a=TRUE, scale_data=TRUE, latent_names=NULL,
                     causal_observed=TRUE, var_observed=TRUE, skew_observed=TRUE, kurt_observed=TRUE,
                     causal_latent=FALSE, var_latent=FALSE, skew_latent=FALSE, kurt_latent=FALSE) {
  # TODO: Expand checks on how many latent can/should be used with or without constrained a depending on input data
  # Input data verification
  if (nrow(data) < 1000)
    stop("Currently only a dataframe with at least 1000 rows is supported.")
  if (ncol(data) < 2) {
    stop("At least two columns in data are required.")
  } else if (ncol(data) == 2) {
    if (n_latent > 1) {
      warning("Only one latent can be used with two variables, unless you really know what you are doing...")
    } else if (!(constrained_a)) {
      warning("Constrained a has to be used with two variables, iunless you really know what you are doing...")
    }
  }
  if (n_latent == 1) {
    if (causal_latent) {stop("latent causal paths are only allowed with >1 latent factor")}
  }
  if (n_latent > ncol(data))
    warning("It's unlikely you want to use use more latent factors than phenotypes present in the data, unless you know what you are doing consider revising....")
  if (any(apply(data, 2, is.character)))
    stop("Numeric column(s) found:",paste0(names(which(apply(data, 2, is.character))), collapse=", "), "\n   Data must only contain numeric columns")
  if (any(round(rowSums(abs(cor(data))), 1) == 1.0)) {
    idcol <- which(round(rowSums(abs(cor(data))), 1) == 1.0)
    stop(paste0("It seems like ",names(idcol), " is an ID column. \n  Consider using ", deparse(substitute(data)), '[,-', unname(idcol),']'))
  }
  if (is.null(latent_names)) {
    if (n_latent > 0) {latent_names <- paste0("f", seq_len(n_latent))} else {latent_names <- as.character(NULL)}
    while (any(latent_names %in% colnames(data))) {
      # Prevent duplicate names between latent and observed
      latent_names[latent_names %in% colnames(data)] <- paste0(latent_names %in% colnames(data), "_0")
    }
  } else {
    if (length(latent_names) != n_latent) {
      stop("Length of latent_names should be equal to n_latent")
    }
  }
  org_names <- if (is.data.frame(data)) {colnames(data)} else {paste0("x", seq_len(ncol(data)))}
  if (is.matrix(data) & !(is.null(colnames(data)))) {org_names <- colnames(data)}
  if (any(org_names %in% latent_names)) {
    stop("Duplicate names found in latent factors and the input dataframe")
  }

  # Scale data
  if (scale_data) {
    data_was_scaled <- all(round(apply(data, 2, mean), 1) == 0) & all(round(apply(data, 2, sd), 1) == 1)
    if (!(data_was_scaled)) {
      data <- apply(data, 2, scale)
    }
  } else {
    data_was_scaled <- TRUE
  }
  data <- as.matrix(data)

  n_p <- ncol(data)
  if (constrained_a) {
    a_names <- NULL
    for (i in seq_len(n_latent)) {
      a_names <- c(a_names, rep(paste0("a", i), n_p))
    }
  } else {
    a_names <- NULL
    for (i in 1:n_latent) {
      for (j in 1:n_p) {
        a_names <- c(a_names, paste0("a", i, "_", j))
      }
    }
  }
  b_names <- NULL
  for (i in 1:n_p) {
    js <- 1:n_p
    for (j in js[js != i]) {
      b_names <- c(b_names, paste0("b",j,"_", i))
    }
  }
  par_names <- list(a=a_names, b=b_names, s=paste0("s", 1:n_p), sk=paste0("sk", 1:n_p), k=paste0("k", 1:n_p))

  named_matrices <- .gen_matrices(par_names, n_p, n_latent, base_value="0")

  ## Make matrices with starting values
  M3.obs <- M3.MM(data)
  M4.obs <- M4.MM(data)
  sk_starts <- NULL
  k_starts <- NULL
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
  num_matrices <- .gen_matrices(par, n_p, n_latent, base_value=0)
  start_values <- as.data.frame(t(unlist(par)))
  colnames(start_values) <- unlist(par_names)
  rownames(start_values) <- "start"
  if ((n_latent > 1) & !(constrained_a)) {
    for (n_fi in 1:(n_latent-1)) {
      num_matrices[["A"]][(n_latent+1):(n_latent+n_fi), n_fi+1] <- 0.0
      named_matrices[["A"]][(n_latent+1):(n_latent+n_fi), n_fi+1] <- "0"
    }
  }
  model <- mcmmodelclass(named_matrices=named_matrices,
                       num_matrices=num_matrices,
                       start_values=start_values,
                       bounds=bounds,
                       meta_data=list(n_obs=nrow(data), n_phenotypes=n_p, n_latent=n_latent, bound_defaults=bound_defaults,
                                      data_was_scaled=data_was_scaled, scale_data=scale_data,
                                      original_colnames=org_names, latent_names=latent_names))
  # causal_observed=TRUE, var_observed=TRUE, skew_observed=TRUE, kurt_observed=TRUE,
  # causal_latent=FALSE, var_latent=FALSE, skew_latent=FALSE, kurt_latent=FALSE
  if (!(causal_observed)) {
    model <- MCMedit(model, "A", 'b', 0)
  }
  if (!(var_observed)) {
    model <- MCMedit(model, "S", 's', 1)
  }
  if (!(skew_observed)) {
    model <- MCMedit(model, "Sk", 'sk', 0)
  }
  if (!(kurt_observed)) {
    model <- MCMedit(model, "K", 'k', 3)
  }

  if (causal_latent) {
    # Add free b-parameters for latent factors
    edit_coords <- list(); edit_coords[[1]] <- edit_coords[[2]] <- -1
    paramnames <- NULL
    for (i in 1:(n_latent-1)) {
      for (j in (i+1):n_latent) {
        edit_coords[[1]] <- c(edit_coords[[1]], i)
        edit_coords[[2]] <- c(edit_coords[[2]], j)
        paramnames <- c(paramnames, paste0("bl",j,"_l",i))
        edit_coords[[1]] <- c(edit_coords[[1]], j)
        edit_coords[[2]] <- c(edit_coords[[2]], i)
        paramnames <- c(paramnames, paste0("bl",i,"_l",j))
      }
    }
    edit_coords[[1]] <- edit_coords[[1]][2:length(edit_coords[[1]])]; edit_coords[[2]] <- edit_coords[[2]][2:length(edit_coords[[2]])]
    model <-  MCMedit(model, "A", edit_coords, paramnames)
  }
  if (var_latent) {
    edit_coords <- list(); edit_coords[[1]] <- edit_coords[[2]] <- -1
    paramnames <- NULL
    for (i in 1:n_latent) {
      edit_coords[[1]] <- c(edit_coords[[1]], i)
      edit_coords[[2]] <- c(edit_coords[[2]], j)
      paramnames <- c(paramnames, paste0("sl",i))
    }
    edit_coords[[1]] <- edit_coords[[1]][2:length(edit_coords[[1]])]; edit_coords[[2]] <- edit_coords[[2]][2:length(edit_coords[[2]])]
    model <- MCMedit(model, "S", edit_coords, paramnames)
  }
  if (skew_latent) {
    edit_coords <- list(); edit_coords[[1]] <- edit_coords[[2]] <- -1
    paramnames <- NULL
    for (i in 1:n_latent) {
      new_coords <- .nd_to_2d_idx(n_latent+n_p, i, i, i)
      edit_coords[[1]] <- c(edit_coords[[1]], new_coords[['x']])
      edit_coords[[2]] <- c(edit_coords[[2]], new_coords[['y']])
      paramnames <- c(paramnames, paste0("skl",i))
    }
    edit_coords[[1]] <- edit_coords[[1]][2:length(edit_coords[[1]])]; edit_coords[[2]] <- edit_coords[[2]][2:length(edit_coords[[2]])]
    model <- MCMedit(model, "Sk", edit_coords, paramnames)
  }
  if (kurt_latent) {
    edit_coords <- list(); edit_coords[[1]] <- edit_coords[[2]] <- -1
    paramnames <- NULL
    for (i in 1:n_latent) {
      new_coords <- .nd_to_2d_idx(n_latent+n_p, i, i, i, i)
      edit_coords[[1]] <- c(edit_coords[[1]], new_coords[['x']])
      edit_coords[[2]] <- c(edit_coords[[2]], new_coords[['y']])
      paramnames <- c(paramnames, paste0("kl",i))
    }
    edit_coords[[1]] <- edit_coords[[1]][2:length(edit_coords[[1]])]; edit_coords[[2]] <- edit_coords[[2]][2:length(edit_coords[[2]])]
    model <- MCMedit(model, "K", edit_coords, paramnames)
  }
  return(model)
}