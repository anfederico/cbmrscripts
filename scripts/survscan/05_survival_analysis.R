km.survfit <- function(x, S, nmin=1, probs=NULL, digits=3, debug=FALSE) {
  if ( !is.Surv(S) )
    stop( "S must be a Surv object" ) 
  if ( nmin>(length(x)-1)/2 )
    stop( paste("nmin too large", nmin) )
  if ( nmin<1 )
    stop( paste("nmin too small:", nmin) )
  if ( !is.null(probs) & nmin>1 )
    warning( "nmin>1 with quantization deprecated" )
  n <- length(unique(x))
  if ( n<nmin*2 )
    stop( paste("not enough x values for given nmin:",n,"(", nmin*2,"needed )") )
  if ( !is.null(probs) & n<length(probs) )
    stop( paste("not enough x values for given quantiles:",n,"(", length(probs),"needed )") )
  
  idx <- !is.na(x)
  x <- x[idx]
  S <- S[idx,]

  # define the set of candidate thresholds
  #
  x.srt <- sort(unique(x))
  x.thresh <- if ( is.null(probs) ) 
    ( x.srt[-1] + x.srt[-length(x.srt)] ) / 2
  else
    quantile(x.srt,probs)
  if (is.null(probs))
    x.thresh <- x.thresh[nmin:(length(x.thresh)-nmin+1)]

  # x.split is a matrix with as many columns as thresholds, with each
  # column indicating which x-values are less than or equal to (FALSE) and
  # greater than (TRUE) the corresponding threshold
  #
  x.split  <- t(apply(matrix(x,length(x),length(x.thresh)),1,'>',x.thresh))

  x.p <- apply( x.split,2,function(z){(1-pchisq(survdiff(S~z)$chisq,1))} )
  best <- which.min( x.p )
  if ( !is.null(digits) ) x.p <- signif(x.p,digits)
  
  if ( !debug )
    return( c( x.thresh[best],      # best trheshold
               x.p[best],           # p-value associated to best threshold
               sum(x.split[,best]), # no. of obs greater than threshold
               sum(S[x.split[,best],1])>sum(S[!x.split[,best],1]))
                                    # does greater than threshold predict good prognosis
           )
  else
    return( cbind(x.thresh, x.p, apply(x.split,2,sum)) )
}

survival.analysis <- function(svdata, time.col, event.col, strata.col, strata.lab, find.thresh=FALSE) {
    oSrv <- Surv(svdata[[time.col]], svdata[[event.col]])
    if (find.thresh) {
        kms.out <- km.survfit(svdata[[strata.col]], S=oSrv, probs=seq(0.05, 0.95, by=0.05), digits=3, debug=TRUE)
        thresh <- kms.out[which.min(kms.out[,'x.p']),'x.thresh']
    } else {
        thresh <- median(svdata[[strata.col]])
    }
    svdata[[strata.lab]] <- ifelse(svdata[[strata.col]] > thresh, 'high', 'low')
    return(list(svdata, thresh))
}

pval.perm <- function(svdata, time.col, event.col, strata.col) {
    X <- svdata[, c(strata.col, strata.col)]
    colnames(X) <- c('c.0', 'c.1')
    X <- t(X)
    oSrv <- Surv(svdata[[time.col]], svdata[[event.col]])
    pout <- perm.survfit(X, surv=oSrv, nperm=1000, probs=seq(0.05,0.95,by=0.05), seed=1234567, verbose=TRUE)
    return(pout)
}

plot_survival_anal <- function(gsva.output, tcga.esets, dir.save=NULL) {

    # Setup dataframes for survival analysis results
    df <- gsva.output
    esets <- tcga.esets
    sigs <- rownames(df[[1]])
    df.srv <- as.data.frame(matrix(0, length(df), length(sigs)))
    rownames(df.srv) <- names(df)
    colnames(df.srv) <- sigs
    df.pval <- df.srv
        
    # Directory for survival plots
    if (!is.null(dir.save)) {
        dir.create(dir.save, recursive=TRUE, showWarnings=FALSE)   
    }
    
    for (i in names(df)) {
        surv_i <- data.frame(t(df[[i]]))
        colnames(surv_i) <- sigs
        pdat_i <- pData(esets[[i]])
        where.subset <- match(rownames(surv_i), rownames(pdat_i))
        surv_i$days_to_death <- pdat_i[where.subset, "days_to_death"]
        surv_i$days_to_last_follow_up <- pdat_i[where.subset, "days_to_last_follow_up"]
        surv_i$death_event <- ifelse(!is.na(surv_i$days_to_death), 1, 0)
        surv_i$days_to_event <- ifelse(!is.na(surv_i$days_to_death), surv_i$days_to_death, surv_i$days_to_last_follow_up)
         
        # Remove NAs
        surv_i <- surv_i[!is.na(surv_i$days_to_event),]      
        time.col <- "days_to_event"
        event.col <- "death_event"
        
        # Store plots
        saved.plots <- list()
    
        for (sig in sigs) {
            strata.col <- sig
            strata.lab <- "group"   
          
            # Determines the stratification threshold
            sv.list <- survival.analysis(surv_i, time.col, event.col, strata.col, strata.lab, find.thresh = TRUE)
            sig.svdata <<- sv.list[[1]] # Global
            sig.svthresh <- sv.list[[2]]
            
            # Survival fit
            svfit <- survfit(Surv(days_to_event, death_event)~group, data=sig.svdata)
            
            # Plotting
            ggsurv <- ggsurvplot(svfit, 
                                pval = TRUE, 
                                xlab = "Time (days)",
                                legend.labs = c(paste("High", svfit$n[1]), paste("Low", svfit$n[2])), 
                                ggtheme = theme_light())
            
            ggsurv$plot <- ggsurv$plot + labs(
              title = "Overall Survival Analysis",                     
              subtitle = paste(i, sig, sep=" : ") 
            )
            
            # Save plot
            saved.plots[[sig]] <- ggsurv$plot
             
            # Save survival results
            surv_i$G <- ifelse(surv_i[,sig] >= mean(surv_i[,sig]), "HI_ACTV", "LO_ACTV")
            pval <- round(surv_pvalue(svfit)$pval, 4)
            survivor <- surv_summary(svfit, surv_i) %>%
                        dplyr::group_by(strata) %>%
                        summarize(msurv=mean(surv))    
            
            survivor.name <- as.character(survivor$strata[which.max(survivor$msurv)])
            df.srv[i, sig] <- survivor.name
            df.pval[i, sig] <- pval
        }
        
        if (!is.null(dir.save)) {
            graphics.off()
            pdf(file.path(dir.save, paste(i, "pdf", sep=".")), onefile=TRUE)
            for (p in saved.plots) {
                show(p)
            }
            graphics.off()        
        }
    }
    
    df.mat <- matrix(0, nrow(df.srv), ncol(df.srv))
    rownames(df.mat) <- rownames(df.srv)
    colnames(df.mat) <- colnames(df.srv)
    
    df.mat[df.srv == "group=high" & df.pval < 0.05] <- 4
    df.mat[df.srv == "group=high" & df.pval >= 0.05] <- 3
    df.mat[df.srv == "group=low" & df.pval >= 0.05] <- 2
    df.mat[df.srv == "group=low" & df.pval < 0.05] <- 1
    
    pvals <- as.numeric(unlist(df.pval[,]))
    ks <- formatC(suppressWarnings(ks.test(pvals, punif, exact=T))$p.value, format="e", digits=2)
    
    order.by <- order(apply(df.mat, 1, sum))
    df.mat.o <- df.mat[order.by,]
    df.pval.o <- df.pval[order.by,]
    
    if (!is.null(dir.save)) {
        graphics.off()
        pdf(file.path(dir.save, paste("Surv-Summary-Plot", "pdf", sep=".")), onefile=TRUE)
        pheatmap(df.mat.o,
                 display_numbers=df.pval.o,
                 cluster_rows=F, 
                 cluster_cols=F,
                 clustering_method="median",
                 labels_col=gsub("\\|", "\n", gsub("_corr", "", colnames(df.mat.o))),
                 fontsize=10,
                 angle_col=0,
                 main=paste0("Survival Analysis\n K.S. P-Value = ", ks),
                 legend=F)
        graphics.off()     
    }
    
    pheatmap(df.mat.o,
             display_numbers=df.pval.o,
             cluster_rows=F, 
             cluster_cols=F,
             clustering_method="median",
             labels_col=gsub("\\|", "\n", gsub("_corr", "", colnames(df.mat.o))),
             fontsize=10,
             angle_col=0,
             main=paste0("Survival Analysis\n K.S. P-Value = ", ks),
             legend=F)
}
