mcmdataclass <- setRefClass("mcmdataclass",
                             fields=list(
                               meta_data="list",
                               M2="matrix",
                               M3="matrix",
                               M4="matrix",
                               SE="list"
                             ))

mcmdataclass$methods(
  initialize=function(meta_data=list(), M2=matrix(), M3=matrix(), M4=matrix(), SE=list()) {
    .self$meta_data <- meta_data
    .self$M2 <- M2
    .self$M3 <- M3
    .self$M4 <- M4
    .self$SE <- SE
  },
  show=function(){
    cat("  MCM data summary\n")
    cat(paste0("N observations : ", .self$meta_data$N, "\n"))
    cat(paste0("N observed vars: ", .self$meta_data$ncol, "\n"))
  },
  copy=function(){
    return(mcmresultclass(.self$meta_data, .self$M2, .self$M3, .self$M4, .self$SE))
  },
  save=function(dest, debug=FALSE) {
    if (!(endsWith(dest, ".mcmdata"))) {dest <- paste0(dest,".mcmdata")}
    file.h5 <- H5File$new(dest, mode="w")
    if (debug) {cat("MCMsavesummary writing M2\n")}
    file.h5$create_dataset("m2", robj = .self$M2, dtype=h5types$H5T_NATIVE_FLOAT, gzip_level = NULL)
    if (debug) {cat("MCMsavesummary writing M3\n")}
    file.h5$create_dataset("m3", robj = .self$M3, dtype=h5types$H5T_NATIVE_FLOAT, gzip_level = NULL)
    if (debug) {cat("MCMsavesummary writing M4\n")}
    file.h5$create_dataset("m4", robj = .self$M4, dtype=h5types$H5T_NATIVE_FLOAT, gzip_level = NULL)
    #file.h5$create_group("data")
    #file.h5[["data/m2"]] <- .self$M2
    #file.h5[["data/m3"]] <- .self$M3
    #file.h5[["data/m4"]] <- .self$M4
    if (debug) {cat("MCMsavesummary writing metadata\n")}
    file.h5$create_group("meta")
    for (i in names(.self$meta_data)) {
      if (is.character(.self$meta_data[[i]])) {
        file.h5[[paste0("meta/", i)]] <- .self$meta_data[[i]]
      } else {
        file.h5[[paste0("meta/", i)]] <- as.numeric(.self$meta_data[[i]])
      }
    }
    file.h5$create_group("SE")
    file.h5[['SE/computed']] <- as.numeric(.self$SE$computed)
    if (.self$SE$computed) {
      if (debug) {cat("MCMsavesummary writing S.m\n")}
      mattotrilvec <- jit_compile(.jit_funcs[['mattotrilvec']])
      #file.h5[['SE/sm']] <- .self$SE$S.m
      file.h5$create_dataset("sm", robj = as.numeric(mattotrilvec$fn(.self$SE$S.m)),
                             dtype=h5types$H5T_NATIVE_FLOAT, gzip_level = NULL)
      if (debug) {cat("MCMsavesummary writing S.m indices\n")}
      if ("idx" %in% names(.self$SE$idx)) {file.h5[['SE/idx']] <- .self$SE$idx$idx}
      if ("idx_nokurt" %in% names(.self$SE$idx)) {file.h5[['SE/idx_nokurt']] <- .self$SE$idx$idx_nokurt}
      if ("idx_noskew" %in% names(.self$SE$idx)) {file.h5[['SE/idx_noskew']] <- .self$SE$idx$idx_noskew}
      file.h5[['SE/idx_nokurt_noskew']] <- .self$SE$idx$idx_nokurt_noskew
    }
    file.h5$close_all()
  },
  load=function(path) {
    file.h5 <- H5File$new(path, mode="r+")
    .self$M2 <- as.matrix(file.h5[["m2"]]$read())
    .self$M3 <- as.matrix(file.h5[["m3"]]$read())
    .self$M4 <- as.matrix(file.h5[["m4"]]$read())
    .self$SE <- list(
      computed=as.logical(file.h5[['SE/computed']]$read()),
      idx=list()
    )
    trilvectomat <- jit_compile(.jit_funcs[['trilvectomat']])
    if (.self$SE$computed) {
      .self$SE[['S.m']] <- trilvectomat$fn(torch_tensor(file.h5[['sm']]$read()))
    }
    idxnames <- list.datasets(file.h5, path="SE/")
    for (idxname in c("idx", "idx_nokurt", "idx_noskew", "idx_nokurt_noskew")) {
      if (idxname %in% idxnames) {
        .self$SE$idx[[idxname]] <- file.h5[[paste0('SE/', idxname)]]$read()
      }
    }
    bools <- c("data_was_scaled", "scale_data", "weighted")  #meta data objects that should be converted to boolean
    for (i in file.h5[["meta"]]$ls()$name) {
      if (i %in% bools) {
        .self$meta_data[[i]] <- as.logical(file.h5[[paste0("meta/", i)]]$read())
      } else {
        .self$meta_data[[i]] <- file.h5[[paste0("meta/", i)]]$read()
      }
    }
    file.h5$close_all()
  }
)
