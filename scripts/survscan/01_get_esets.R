get_esets <- function(tumor.types=c(1), # Primary Solid Tumor
                      normalization.method="DESeq2", # c("DESeq2", "edgeR")
                      dir.cache=NULL) {

    pattern <- gsub("\\{0}", normalization.method, "*{0}_log_eSet.rds")
    files <- list.files(path=file.path(path.to.tcgadump, "esets_normalized"),
                        pattern=pattern, full.names=TRUE, recursive=FALSE)
    
    esets <- list()
    for (fn in files) {
          splitted <- strsplit(fn, "/")[[1]]
          id <- strsplit(splitted[[length(splitted)]], "_")[[1]][1]
          eset <- readRDS(fn)
          
          # We just want Primary Solid Tumor (Code = 1)
          eset <- eset[,eset$sample_type_id %in% tumor.types]
         
          # Make sure survival information is correctly formatted
          pdat <- pData(eset)
          
          # Suppress warnings because its introducing NAs to missing entries
          suppressWarnings(pdat$days_to_death <- as.numeric(pdat$days_to_death))
          suppressWarnings(pdat$days_to_last_follow_up <- as.numeric(pdat$days_to_last_follow_up))
          
          pData(eset) <- pdat
          esets[[id]] <- eset
        
          rm(splitted, id, eset, pdat) # Cleanup
    }
    rm(pattern, fn, files)
    
    # Remove datasets without solid primary tumor samples
    DIM <- data.frame(t(sapply(esets,dim)))
    invalid <- rownames(DIM)[DIM$Samples < 1]
    for (i in invalid) {
        cat("Removing...", i, "\n")
        esets[[i]] <- NULL
    }
    names(esets)
    if (!is.null(dir.cache)) {
        dir.create(file.path(dir.cache, "cache"), showWarnings=FALSE)
        saveRDS(esets, file.path(dir.cache, "cache/tcga.esets.rds"))
    }
    
    DIM <- data.frame(t(sapply(esets,dim)))
    DIM$`Has Days to Death` <- as.numeric(lapply(esets, function(x) table(is.na(pData(x)$days_to_death))["FALSE"]))
    DIM$`Has Days to Last Followup` <- as.numeric(lapply(esets, function(x) table(is.na(pData(x)$days_to_last_follow_up))["FALSE"]))
    rownames(DIM) <- names(esets)
    print(DIM)  
    
    return(esets)
}
