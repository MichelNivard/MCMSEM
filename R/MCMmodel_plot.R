plot.mcmmodelclass <- function(model,layout = NULL, use_values=FALSE, ...) {
  #TODO: I think this is not required, remove it in the next version if nothing weird happens
  #latents <- seq_len(model$meta_data$n_latent)
  #if (model$meta_data$n_latent >= 1) {
  #  qgraph(t(model$num_matrices$A[-latents,-latents] + model$num_matrices$S[-latents,-latents]), ...)
  #} else {
  #  qgraph(t(model$num_matrices$A + model$num_matrices$S), ...)
  #}
  # Full model
  place.latents <- floor(quantile(seq_len(model$meta_data$n_phenotypes),seq_len(model$meta_data$n_latent)/(model$meta_data$n_latent+1)))
  if (model$meta_data$n_latent == 1) place.latents <- place.latents + 0.5
  if (model$meta_data$n_phenotypes < 10) {
    pheno_pos <-rep(1,model$meta_data$n_phenotypes)
  } else {
    pheno_pos <- NULL
    for (i in 1:model$meta_data$n_phenotypes) {
      pheno_pos <- c(pheno_pos, ((i %% 2) / 2) + ((i-1) %% 2))  # simplest way I knew how to get 1.0, 0.5, 1.0, 0.5, ....
    }
  }
  if(is.null(layout)){
  layout <- cbind(c(place.latents,1:model$meta_data$n_phenotypes),
                  c(rep(2,model$meta_data$n_latent),pheno_pos))
    }
  plot_mat <- model$num_matrices$A[,] + model$num_matrices$S[,]
  if (!(use_values)) {
    which_0 <- (model$named_matrices$A == "0") & (model$named_matrices$S == "0")
    plot_mat[,] <- 1.0
    # If all the variances are equal they are not plotted, so I make them every so slightly different to hack them into the qgraph
    diag(plot_mat) <- rnorm(nrow(model$named_matrices$A ), 1, sd=.Machine$double.eps)
    plot_mat[which_0] <- 0.0
  }
  rownames(plot_mat) <- c(model$meta_data$latent_names, model$meta_data$original_colnames)
  colnames(plot_mat) <- c(model$meta_data$latent_names, model$meta_data$original_colnames)
  if (use_values) {
    return(qgraph(t(plot_mat),layout=layout, ...))
  } else {
    return(qgraph(t(plot_mat),layout=layout, minimum=0, maximum=1.01, mode="direct", weighted=FALSE, ...))
  }
}
