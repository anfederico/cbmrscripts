combine_signatures <- function(gsva.output, pairs, keep.original=FALSE) {
    lapply(gsva.output, function(x) {
        for (p in pairs) {
            row <- x[p[2],]-x[p[3],]
            x <- rbind(x, row)
            rownames(x) <- c(rownames(x)[1:nrow(x)-1], p[1])
        }
        if (keep.original) {
            return(x)
        } 
        else { # Drop old signatures
            combined <- unlist(lapply(pairs, function(x) x[1]))
            return(x[combined,,drop=FALSE])
        }
    })
}
