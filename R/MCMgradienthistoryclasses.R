## MCM gradienthistory class
mcmgradienthistoryclass <- setRefClass("mcmgradienthistoryclass",
                             fields=list(
                               df="data.frame",
                               label="character",
                               hasgrads="logical",
                               last_iter="data.frame"
                             ))

mcmgradienthistoryclass$methods(
  initialize=function(df=data.frame(NULL), label=as.character(NULL), hasgrads=FALSE, last_iter=data.frame(NULL)) {
    .self$df <- df
    .self$label <- label
    .self$hasgrads <- hasgrads
    .self$last_iter <- last_iter
  },
  show=function(){
    if (.self$hasgrads) {
      cat(paste0("  Gradient history ",.self$label," parameters\n"))
      print(summary(.self))
    } else {
      cat(paste0("No parameters that require gradients in ", .self$label, "\n"))
    }
  },
  copy=function(){
    return(mcmgradienthistoryclass(.self$df, .self$label, .self$hasgrads, .self$last_iter))
  }
)

as.data.frame.mcmgradienthistoryclass <- function(x) {
  return(x$df)
}

as.matrix.mcmgradienthistoryclass <- function(x) {
  return(as.matrix(x$df))
}

abs.mcmgradienthistoryclass <- function(x) {
  return(mcmgradienthistoryclass(abs(x$df), label=x$label, hasgrads=x$hasgrads, last_iter=abs(x$last_iter)))
}

log.mcmgradienthistoryclass <- function(x, base=exp(1)) {
  if (any(x$df < 0)) {
    stop("Error: negative values in MCM gradient history. Consider using log(abs(x))")
  }
  if (any(x$df == 0)) {
    warning("Values of 0 detected, adding minimum ammount to all values to prevent infinites")
    return(mcmgradienthistoryclass(log(x$df+.Machine$double.eps, base=base), label=x$label, hasgrads=x$hasgrads, last_iter=data.frame(log(x$last_iter, base=base))))
  }
  return(mcmgradienthistoryclass(log(x$df, base=base), label=x$label, hasgrads=x$hasgrads, last_iter=data.frame(log(x$last_iter, base=base))))
}

log10.mcmgradienthistoryclass <- function(x) {
  return(log(x, base=10))
}

sqrt.mcmgradienthistoryclass <- function(x) {
  if (any(x$df < 0)) {
    stop("Error: negative values in MCM gradient history. Consider using sqrt(abs(x))")
  }
  return(mcmgradienthistoryclass(sqrt(x$df), label=x$label, hasgrads=x$hasgrads, last_iter=data.frame(sqrt(x$last_iter))))
}

min.mcmgradienthistoryclass <- function(x, na.rm=TRUE) {
  return(apply(x$df, 2, min, na.rm=na.rm))
}

max.mcmgradienthistoryclass <- function(x, na.rm=TRUE) {
  return(apply(x$df, 2, max, na.rm=na.rm))
}

summary.mcmgradienthistoryclass <- function(x) {
  if (x$hasgrads) {
    y <- rbind(min(x), max(x), min(abs(x)), max(abs(x)), x$last_iter)
    rownames(y) <- c("min", "max", "min(abs)", "max(abs)", "last iteration")
    return(y)
  } else {
    return(data.frame(NULL))
  }
}

plot.mcmgradienthistoryclass <- function(x, parameterlegend=TRUE, xlim='auto', ylim='auto', cols=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E142", "#0072B2", "#D55E00", "#CC79A7"),
                                         ltys=c("solid", 'dashed', 'dotted', 'dotdash'), main='Gradient history of %s parameters',xlab='iteration', ylab='gradient', legendloc='topright', bty='n', ...) {
  # by default: 8 cols * 4 ltys to allow for 32 unique combinations
  if (x$hasgrads) {
    if (all(xlim == 'auto')) {xlim <- c(0, nrow(x$df))}
    if (all(ylim == 'auto')) {
      if (min(x$df) >= 0) {
        ylim <- c(0, max(abs(x$df)))
      } else {
        ylim <- c(-max(abs(x$df)), max(abs(x$df)))
      }
    }
    if (grepl('%s', main)) {main <- sprintf(main, x$label)}
    plot(-100, -100, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, bty=bty, ...)
    legend_txt <- c()
    legend_col <- c()
    legend_lty <- c()
    for (npar in seq_len(ncol(x$df))) {
      colname <- colnames(x$df)[npar]
      col <- cols[((npar-1) %% length(cols)) + 1]
      lty <- ltys[(floor((npar-1)/length(ltys)) %% length(ltys)) + 1]
      lines(seq_len(nrow(x$df)), x$df[, npar],
            col=col,
            lty=lty)
      legend_txt <- c(legend_txt, colname)
      legend_col <- c(legend_col, col)
      legend_lty <- c(legend_lty, lty)
    }
    if (parameterlegend) {
      legend(legendloc, legend=legend_txt, col=legend_col, lty=legend_lty)
    }
  } else {
    plot(-100, -100, xlim=c(0, 1), ylim=c(0, 1), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
    text(0.5, 0.5, paste0("No parameters that require gradients in ", .self$label))
  }
}

## MCM multi-gradienthistory class
mcmmultigradienthistoryclass <- setRefClass("mcmmultigradienthistoryclass",
                             fields=list(
                               A="mcmgradienthistoryclass",
                               Fm="mcmgradienthistoryclass",
                               S="mcmgradienthistoryclass",
                               Sk="mcmgradienthistoryclass",
                               K="mcmgradienthistoryclass",
                               hasgrads="logical"
                             ))

mcmmultigradienthistoryclass$methods(
  initialize=function(x=list(), A=mcmgradienthistoryclass(label='A'), Fm=mcmgradienthistoryclass(label='Fm'), S=mcmgradienthistoryclass(label='S'),
                      Sk=mcmgradienthistoryclass(label='Sk'), K=mcmgradienthistoryclass(label='K'), hasgrads=FALSE) {
    if (length(x) > 0) {
      .self$hasgrads <- hasgrads
      for (matname in c("A", "Fm", "S", "Sk", "K")) {
        if (matname %in% names(x)) {
          if (nrow(x[[matname]]) > 1) {
            .self[[matname]] <- mcmgradienthistoryclass(x[[matname]], label=matname, hasgrads=TRUE, last_iter=x[[matname]][nrow(x[[matname]]), ])
          } else {
            .self[[matname]] <- mcmgradienthistoryclass(data.frame(NULL), label=matname, hasgrads=FALSE, last_iter=x[[matname]][1, ])
          }
        } else {
          .self[[matname]] <- mcmgradienthistoryclass(label=matname, hasgrads=FALSE)
        }
      }
    } else {
      .self$A  <- A
      .self$Fm <- Fm
      .self$S  <- S
      .self$Sk <- Sk
      .self$K  <- K
      .self$hasgrads <- hasgrads
    }
  },
  show=function(){
    if (.self$hasgrads) {
      cat(paste0("  Gradient history summary\n"))
      for (i in c("A", "Fm", "S", "Sk", "K")) {
        if (.self[[i]]$hasgrads) {
          print(summary(.self[[i]]))
        }
      }
    } else {
      cat("Gradient history not stored, run MCMfit with monitor_grads=TRUE\n")
    }
  },
  copy=function(){
    return(A=.self$A, Fm=.self$Fm, S=.self$S, Sk=.self$Sk, K=.self$K, hasgrads=.self$hasgrads)
  }
)

as.data.frame.mcmmultigradienthistoryclass <- function(x) {
  if (x$hasgrads) {
    res <- list()
    for (i in c("A", "Fm", "S", "Sk", "K")) {
      if (x[[i]]$hasgrads) {
        res[[i]] <- x[[i]]$df
      }
    }
    return(do.call(cbind, res))
  } else {
    return(data.frame(NULL))
  }
}

as.matrix.mcmmultigradienthistoryclass <- function(x) {
  return(as.matrix(as.data.frame(x)))
}

abs.mcmmultigradienthistoryclass <- function(x) {
  return(mcmmultigradienthistoryclass(
    A=abs(x$A), Fm=abs(x$Fm), S=abs(x$S), Sk=abs(x$Sk), K=abs(x$K), hasgrads=x$hasgrads
  ))
}

log.mcmmultigradienthistoryclass <- function(x, base=exp(1)) {
  return(mcmmultigradienthistoryclass(
    A=log(x$A, base=base), Fm=log(x$Fm, base=base), S=log(x$S, base=base), Sk=log(x$Sk,base=base), K=log(x$K,base=base), hasgrads=x$hasgrads
  ))
}

log10.mcmmultigradienthistoryclass <- function(x, base=exp(1)) {
  return(log(x, base=10))
}

sqrt.mcmmultigradienthistoryclass <- function(x) {
  return(mcmmultigradienthistoryclass(
    A=sqrt(x$A), Fm=sqrt(x$Fm), S=sqrt(x$S), Sk=sqrt(x$Sk), K=sqrt(x$K), hasgrads=x$hasgrads
  ))
}

min.mcmmultigradienthistoryclass <- function(x, na.rm=TRUE) {
  if (x$hasgrads) {
    res <- c()
    for (i in c("A", "Fm", "S", "Sk", "K")) {
      if (x[[i]]$hasgrads) {
        res <- c(res, min(x[[i]], na.rm=na.rm))
      }
    }
    return(res)
  } else {
    return(NULL)
  }
}

max.mcmmultigradienthistoryclass <- function(x, na.rm=TRUE) {
  if (x$hasgrads) {
    res <- c()
    for (i in c("A", "Fm", "S", "Sk", "K")) {
      if (x[[i]]$hasgrads) {
        res <- c(res, max(x[[i]], na.rm=na.rm))
      }
    }
    return(res)
  } else {
    return(NULL)
  }
}

summary.mcmmultigradienthistoryclass <- function(x) {
  if (x$hasgrads) {
    res <- list()
    for (i in c("A", "Fm", "S", "Sk", "K")) {
      if (x[[i]]$hasgrads) {
        res[[i]] <- summary(x[[i]])
      }
    }
    return(do.call(cbind, res))
  } else {
    return(data.frame(NULL))
  }
}

plot.mcmmultigradienthistoryclass <- function(x, layout.matrix='auto', parameterlegend=TRUE, xlim='auto', ylim='auto', cols=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E142", "#0072B2", "#D55E00", "#CC79A7"),
                                         ltys=c("solid", 'dashed', 'dotted', 'dotdash'), main='Gradient history of %s parameters',xlab='iteration', ylab='gradient', legendloc='topright', bty='n', ...) {
  if (x$hasgrads) {
    pars_to_plot <- c()
    for (i in c("A", "Fm", "S", "Sk", "K")) {
      if (x[[i]]$hasgrads) {
        pars_to_plot <- c(pars_to_plot, i)
      }
    }
    if (all(layout.matrix == 'auto')) {
      if (length(pars_to_plot) == 5) {layout.matrix <- matrix(c(1, 1, 2:5), ncol=2, byrow=TRUE)}
      if (length(pars_to_plot) == 4) {layout.matrix <- matrix(1:4, ncol=2, byrow=TRUE)}
      if (length(pars_to_plot) == 3) {layout.matrix <- matrix(c(1, 1, 2:3), ncol=2, byrow=TRUE)}
      if (length(pars_to_plot) == 2) {layout.matrix <- matrix(c(1, 2), ncol=1, byrow=TRUE)}
      if (length(pars_to_plot) == 1) {layout.matrix <- matrix(1, ncol=1, byrow=TRUE)}
    } else {
      if (length(unique(layout.matrix)) != length(pars_to_plot)) {
        stop(paste0("Expected a layout matrix with ", length(pars_to_plot)," unique values"))
      }
    }
    layout(layout.matrix)
    for (i in pars_to_plot) {
      plot(x[[i]], parameterlegend=parameterlegend, xlim=xlim, ylim=ylim, cols=cols, ltys=ltys, main=main, xlab=xlab, ylab=ylab, legendloc=legendloc, bty=bty)
    }
  } else {
    plot(-100, -100, xlim=c(0, 1), ylim=c(0, 1), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
    text(0.5, 0.5, "Gradient history not stored, run MCMfit with monitor_grads=TRUE")
  }
}