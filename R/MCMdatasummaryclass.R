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
  save=function(dest) {
    if (!(endsWith(dest, ".mcmdata"))) {dest <- paste0(dest,".mcmdata")}
    file.h5 <- H5File$new(dest, mode="w")
    file.h5$create_dataset("m2", robj = .self$M2, chunk_dims = "auto", gzip_level = 9)
    file.h5$create_dataset("m3", robj = .self$M3, chunk_dims = "auto", gzip_level = 9)
    file.h5$create_dataset("m4", robj = .self$M4, chunk_dims = "auto", gzip_level = 9)
    #file.h5$create_group("data")
    #file.h5[["data/m2"]] <- .self$M2
    #file.h5[["data/m3"]] <- .self$M3
    #file.h5[["data/m4"]] <- .self$M4
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
      #file.h5[['SE/sm']] <- .self$SE$S.m
      file.h5$create_dataset("sm", robj = .self$SE$S.m, chunk_dims = "auto", gzip_level = 9)
      file.h5[['SE/idx']] <- .self$SE$idx$idx
      file.h5[['SE/idx_nokurt']] <- .self$SE$idx$idx_nokurt
      file.h5[['SE/idx_noskew']] <- .self$SE$idx$idx_noskew
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
      S.m=file.h5[['sm']]$read(),
      idx=list(
        idx=file.h5[['SE/idx']]$read(),
        idx_nokurt=file.h5[['SE/idx_nokurt']]$read(),
        idx_noskew=file.h5[['SE/idx_noskew']]$read(),
        idx_nokurt_noskew=file.h5[['SE/idx_nokurt_noskew']]$read()
      )
    )
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
