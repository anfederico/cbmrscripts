do_gsva <- function(signatures, tcga.esets, dir.cache=NULL) {
    
    datasets <- names(tcga.esets)
    gsva.output <- list()
    for (i in datasets) {
      cat(i, "...\n")
      gsva.output[[i]] <- GSVA::gsva(exprs(tcga.esets[[i]]), signatures, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
    }

    if (!is.null(dir.cache)) {
        dir.create(file.path(dir.cache, "cache"), showWarnings = FALSE)
        saveRDS(gsva.output, file.path(dir.cache, "cache/gsva.output.rds"))
    }

    return(gsva.output)
}