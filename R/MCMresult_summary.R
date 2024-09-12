summary.mcmresultclass <- function(res) {
  Pars_reg <- data.frame(matrix(NA, ncol=8, nrow=1))
  Pars_fact <- data.frame(matrix(NA, ncol=8, nrow=1))
  for (col in seq_len(ncol(res$model$named_matrices[['A']]))) {
    for (row in seq_len(nrow(res$model$named_matrices[['A']]))) {
      if (is.na(suppressWarnings(as.numeric(res$model$named_matrices[['A']][row, col])))) {
        parname <- res$model$named_matrices[['A']][row, col]
        parvalue <- res$model$num_matrices[['A']][row, col]
        std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
        p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
        if (xor(row > res$model$meta_data$n_latent, col > res$model$meta_data$n_latent)) {
          if (col <= res$model$meta_data$n_latent) {
            lhs <- res$model$meta_data$latent_names[col]
            rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_latent]
            Pars_reg <- rbind(Pars_reg, c(parname, lhs, "=~", rhs, parvalue, std, p, as.numeric(res$gradients$A$last_iter[parname])))
          } else {
            lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_latent]
            rhs <- res$model$meta_data$latent_names[row]
            Pars_reg <- rbind(Pars_reg, c(parname, lhs, "~>", rhs, parvalue, std, p, as.numeric(res$gradients$A$last_iter[parname])))
          }
        } else {
          if (row > res$model$meta_data$n_latent) {
            lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_latent]
            rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_latent]
          } else {
            lhs <- res$model$meta_data$latent_names[col]
            rhs <- res$model$meta_data$latent_names[row]
          }
          Pars_reg <- rbind(Pars_reg, c(parname, lhs, "~>", rhs, parvalue, std, p, as.numeric(res$gradients$A$last_iter[parname])))
        }
      }
    }
  }
  # Drop the first row safely, i.e. if there is only the NA row, this results in
  Pars_fact <- Pars_fact[seq_len(nrow(Pars_fact)-1)+1,]
  Pars_reg <- Pars_reg[seq_len(nrow(Pars_reg)-1)+1,]
  Pars <- rbind(Pars_fact, Pars_reg)
  colnames(Pars) <- c("label", "lhs", "edge", "rhs", "est", "se", "p", "last_gradient")
  rownames(Pars) <- seq_len(nrow(Pars))
  Vars <- data.frame(matrix(NA, ncol=8, nrow=1))
  for (row in seq_len(nrow(res$model$named_matrices[['S']]))) {
    for (col in seq_len(ncol(res$model$named_matrices[['S']]))) {
      if (is.na(suppressWarnings(as.numeric(res$model$named_matrices[['S']][row, col])))) {
        parname <- res$model$named_matrices[['S']][row, col]
        parvalue <- res$model$num_matrices[['S']][row, col]
        if (col <= res$model$meta_data$n_latent) {
          lhs <- res$model$meta_data$latent_names[col]
        } else {
          lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_latent]
        }
        if (row <= res$model$meta_data$n_latent) {
          rhs <- res$model$meta_data$latent_names[row]
        } else {
          rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_latent]
        }
        label <- parname
        edge <- "~~"
        est <- parvalue
        std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
        p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
        Vars <- rbind(Vars, c(label, lhs, edge, rhs, est, std, p, as.numeric(res$gradients$S$last_iter[parname])))
      }
    }
  }
  Vars <- Vars[2:nrow(Vars), ]
  colnames(Vars) <- c("label", "lhs", "edge", "rhs", "est", "se", "p", "last_gradient")
  rownames(Vars) <- seq_len(nrow(Vars))
  if (res$info$use_skewness) {
    Skews <- data.frame(matrix(NA, ncol=9, nrow=1))
    for (row in seq_len(nrow(res$model$named_matrices[['Sk']]))) {
      for (col in seq_len(ncol(res$model$named_matrices[['Sk']]))) {
        if (is.na(suppressWarnings(as.numeric(res$model$named_matrices[['Sk']][row, col])))) {
          parname <- res$model$named_matrices[['Sk']][row, col]
          parvalue <- res$model$num_matrices[['Sk']][row, col]
          vars <- .twod_to_nd_idx(nrow(res$model$named_matrices[['Sk']]), row, col, ndims=3)
          v1 <- if (vars[[1]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[1]]]} else {res$model$meta_data$original_colnames[vars[[1]] - res$model$meta_data$n_latent]}
          v2 <- if (vars[[2]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[2]]]} else {res$model$meta_data$original_colnames[vars[[2]] - res$model$meta_data$n_latent]}
          v3 <- if (vars[[3]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[3]]]} else {res$model$meta_data$original_colnames[vars[[3]] - res$model$meta_data$n_latent]}
          label <- parname
          edge <- "~~~"
          est <- parvalue
          std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
          p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
          Skews <- rbind(Skews, c(label, edge, v1, v2, v3, est, std, p, as.numeric(res$gradients$Sk$last_iter[parname])))
        }
      }
    }
    Skews <- Skews[2:nrow(Skews), ]
    colnames(Skews) <- c("label", "edge", "v1", "v2", "v3", "est", "se", "p", "last_gradient")
    rownames(Skews) <- seq_len(nrow(Skews))
  }
  if (res$info$use_kurtosis) {
    Kurts <- data.frame(matrix(NA, ncol=10, nrow=1))
    for (row in seq_len(nrow(res$model$named_matrices[['K']]))) {
      for (col in seq_len(ncol(res$model$named_matrices[['K']]))) {
        if (is.na(suppressWarnings(as.numeric(res$model$named_matrices[['K']][row, col])))) {
          parname <- res$model$named_matrices[['K']][row, col]
          parvalue <- res$model$num_matrices[['K']][row, col]
          vars <- .twod_to_nd_idx(nrow(res$model$named_matrices[['K']]), row, col, ndims=4)
          v1 <- if (vars[[1]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[1]]]} else {res$model$meta_data$original_colnames[vars[[1]] - res$model$meta_data$n_latent]}
          v2 <- if (vars[[2]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[2]]]} else {res$model$meta_data$original_colnames[vars[[2]] - res$model$meta_data$n_latent]}
          v3 <- if (vars[[3]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[3]]]} else {res$model$meta_data$original_colnames[vars[[3]] - res$model$meta_data$n_latent]}
          v4 <- if (vars[[4]] <= res$model$meta_data$n_latent) {res$model$meta_data$latent_names[vars[[4]]]} else {res$model$meta_data$original_colnames[vars[[4]] - res$model$meta_data$n_latent]}
          label <- parname
          edge <- "~~~~"
          est <- parvalue
          std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
          p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
          Kurts <- rbind(Kurts, c(label, edge, v1, v2, v3, v4, est, std, p, as.numeric(res$gradients$K$last_iter[parname])))
        }
      }
    }
    Kurts <- Kurts[2:nrow(Kurts), ]
    colnames(Kurts) <- c("label", "edge", "v1", "v2", "v3", "v4", "est", "se", "p", "last_gradient")
    rownames(Kurts) <- seq_len(nrow(Kurts))
  }
  for (i in c("est", "se", "p")) {
    Pars[, i] <- as.numeric(Pars[, i])
    Vars[, i] <- as.numeric(Vars[, i])
    if (res$info$use_skewness) {Skews[, i] <- as.numeric(Skews[, i])}
    if (res$info$use_kurtosis) {Kurts[, i] <- as.numeric(Kurts[, i])}
  }

  loss <- res$loss
  n_par <- length(res$model$param_values)
  n_obs <- res$model$meta_data$n_obs
  chisq <- n_obs * loss
  # loosely according to Broudt et al:
  # AIC disabled for now aic <- n_obs * loss + 2 * n_par
  bic <- n_obs * loss + n_par * log(n_obs)
  return(mcmresultsummaryclass(parameters=Pars, variances=Vars,
                               skewness=if (res$info$use_skewness) {Skews} else {as.data.frame(NULL)},
                               kurtosis=if (res$info$use_kurtosis) {Kurts} else {as.data.frame(NULL)},
                               loss=loss, n_par=n_par, n_obs=n_obs, chisq=chisq, bic=bic, result=res$copy()))
}