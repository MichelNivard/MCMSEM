summary.mcmresultclass <- function(res) {
  Pars_reg <- data.frame(matrix(NA, ncol=9, nrow=1))
  Pars_fact <- data.frame(matrix(NA, ncol=9, nrow=1))
  iter <- 1
  for (col in 1:ncol(res$model$named_matrices[['A']])) {
    for (row in 1:nrow(res$model$named_matrices[['A']])) {
      if (res$model$named_matrices[['A']][row, col] != "0") {
        parname <- res$model$named_matrices[['A']][row, col]
        parvalue <- res$model$num_matrices[['A']][row, col]
        if (startsWith(parname, "a")) {
          label <- parname
          if (col < row) {
            lhs <- res$model$meta_data$latent_names[col]
            rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_confounding]
          } else {
            lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_confounding]
            rhs <- res$model$meta_data$latent_names[row]
          }
          edge <- "=~"
          est <- parvalue
          std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
          p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
          group <- 1; fixed <- FALSE
          Pars_fact <- rbind(Pars_fact, c(label, lhs, edge, rhs, est, std, p, group, fixed))
        } else {
          label <- parname
          lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_confounding]
          rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_confounding]
          edge <- "~>"
          est <- parvalue
          std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
          p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
          group <- 1; fixed <- FALSE
          Pars_reg <- rbind(Pars_reg, c(label, lhs, edge, rhs, est, std, p, group, fixed))
        }
        iter <- iter + 1
      }
    }
  }
  Pars_fact <- Pars_fact[2:nrow(Pars_fact), ]
  Pars_reg <- Pars_reg[2:nrow(Pars_reg), ]
  Pars <- rbind(Pars_fact, Pars_reg)
  for (row in 1:nrow(res$model$named_matrices[['S']])) {
    for (col in 1:ncol(res$model$named_matrices[['S']])) {
      if (!(res$model$named_matrices[['S']][row, col] %in% c("0", "1"))) {
        parname <- res$model$named_matrices[['S']][row, col]
        parvalue <- res$model$num_matrices[['S']][row, col]
        if (col <= res$model$meta_data$n_confounding) {
          lhs <- res$model$meta_data$latent_names[col]
        } else {
          lhs <- res$model$meta_data$original_colnames[col - res$model$meta_data$n_confounding]
        }
        if (row <= res$model$meta_data$n_confounding) {
          rhs <- res$model$meta_data$latent_names[row]
        } else {
          rhs <- res$model$meta_data$original_colnames[row - res$model$meta_data$n_confounding]
        }
        label <- parname
        edge <- "~~"
        est <- parvalue
        std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
        p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
        group <- 1; fixed <- FALSE
        Pars <- rbind(Pars, c(label, lhs, edge, rhs, est, std, p, group, fixed))
        iter <- iter + 1
      }
    }
  }
  # TODO: For now Sk/K output is fixed, instead of obtaining the parameter names/values from the actual matrix...
  #  in short, it will work with default settings, but if users add new co-kurtosis/co-skewness parameters these will not be included
  for (i in 1:res$model$meta_data$n_phenotypes) {
    parname <- paste0("sk", i)
    if (parname %in% colnames(res$df)) {
      parvalue <- res$df['est', parname]
      lhs <- res$model$meta_data$original_colnames[i]
      rhs <- res$model$meta_data$original_colnames[i]
      label <- parname
      edge <- "~~~"
      est <- parvalue
      std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
      p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
      group <- 1; fixed <- FALSE
      Pars <- rbind(Pars, c(label, lhs, edge, rhs, est, std, p, group, fixed))
      iter <- iter + 1
    }
  }
  for (i in 1:res$model$meta_data$n_phenotypes) {
    parname <- paste0("k", i)
    if (parname %in% colnames(res$df)) {
      parvalue <- res$df['est', parname]
      lhs <- res$model$meta_data$original_colnames[i]
      rhs <- res$model$meta_data$original_colnames[i]
      label <- parname
      edge <- "~~~~"
      est <- parvalue
      std <- if ('se' %in% rownames(res$df)) {res$df['se', parname]} else {NA}
      p <- if ('se' %in% rownames(res$df)) {2*pnorm(abs(res$df['est', parname])/res$df['se', parname], lower.tail=FALSE)} else {NA}
      group <- 1; fixed <- FALSE
      Pars <- rbind(Pars, c(label, lhs, edge, rhs, est, std, p, group, fixed))
      iter <- iter + 1
    }
  }

  colnames(Pars) <- c("label", "lhs", "edge", "rhs", "est", "se", "p", "group", "fixed")
  Pars$est <- as.numeric(Pars$est)
  Pars$se <- as.numeric(Pars$se)
  Pars$p <- as.numeric(Pars$p)
  rownames(Pars) <- 1:nrow(Pars)
  loss <- res$loss
  n_par <- length(res$model$param_values)
  n_obs <- res$model$meta_data$n_obs
  chisq <- n_obs * loss
  # loosely according to Broudt et al:
  aic <- n_obs * loss + 2 * n_par
  bic <- n_obs * loss + n_par * log(n_obs)
  return(mcmresultsummaryclass(df=Pars, loss=loss, n_par=n_par, n_obs=n_obs, chisq=chisq, aic=aic, bic=bic, result=res))
}