add_timer <- function(gsva.output, tcga.esets) {
    datasets <- names(tcga.esets)
    timer <- read.table(path.to.timer, header=T, sep="\t")
    df <- gsva.output
    for (i in datasets) {
        df_i <- as.data.frame(t(df[[i]]))
        sample.ids <- substr(rownames(df_i), 1, 16)
        where.subset <- match(sample.ids, timer$SampleID)
        df_i$CD4_INFL <- timer[where.subset, 'T_cell.CD4']
        df_i$CD8_INFL <- timer[where.subset, 'T_cell.CD8']
        df_i$SUM_INFL <- df_i$CD4_INFL+df_i$CD8_INFL
        df[[i]] <- df_i
    }
    return(df)  
}

plot_immune_corr <- function(gsva.output, tcga.esets, dir.save=NULL) {
    sig_names <- rownames(gsva.output[[1]])
    df <- add_timer(gsva.output, tcga.esets)
    
    # Setup dataframe
    df.infl <- as.data.frame(matrix(0, nrow=length(df), ncol=length(sig_names)*3+1))
    cn <- sapply(sig_names, function(x) paste(x, 'corr', sep="_"))
    pn <- sapply(sig_names, function(x) paste(x, 'pval', sep="_"))
    rn <- sapply(sig_names, function(x) paste(x, 'r2', sep="_"))
    colnames(df.infl) <- c(c('SUM_INFL'), as.character(rbind(cn, pn, rn)))
    rownames(df.infl) <- names(df)
    for (i in rownames(df.infl)) {
        df.infl[i, 'SUM_INFL'] <- mean(df[[i]]$SUM_INFL, na.rm=TRUE)
    }
    
    # Remove TCGA datasets without infiltration data
    df.infl <- na.omit(df.infl)
      
    # Correlation of signatures with CD4 and CD8 Infiltration
    for (i in rownames(df.infl)) {
        df_i <- df[[i]]
        for (sig in sig_names) {
            fit = lm(df_i[['SUM_INFL']] ~ df_i[[sig]], df_i)
            df.infl[i, paste(sig, "corr", sep="_")] <- summary(fit)$coefficients[2,1]
            df.infl[i, paste(sig, "pval", sep="_")] <- summary(fit)$coefficients[2,4]  
            df.infl[i, paste(sig, "r2", sep="_")] <- summary(fit)$r.squared
        }
    }
    
    df.corr <- df.infl[,seq(2, ncol(df.infl), by=3)][]
    df.pval <- df.infl[,seq(3, ncol(df.infl), by=3)]
    df.pval[] <- lapply(df.pval, function(x) formatC(x, format='e', 2))
    
    if (!is.null(dir.save)) {
        graphics.off()
        pdf(dir.save, onefile=TRUE)
        pheatmap(df.corr,
                 display_numbers=df.pval,
                 cluster_rows=T, 
                 cluster_cols=F,
                 fontsize_number=10, 
                 labels_col=gsub("\\|", "\n", gsub("_corr", "", colnames(df.corr))),
                 fontsize=10,
                 angle_col=0,
                 main="Correlation of Signatures with Infiltration")
        graphics.off()
    }

    # Show in session
    pheatmap(df.corr,
             display_numbers=df.pval,
             cluster_rows=T, 
             cluster_cols=F,
             fontsize_number=10, 
             labels_col=gsub("\\|", "\n", gsub("_corr", "", colnames(df.corr))),
             fontsize=10,
             angle_col=0,
             main="Correlation of Signatures with Infiltration")
}
