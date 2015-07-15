filterRes <- function(res,
                      gp1, gp2, 
                      background_coverage_threshold=2, 
                      P.Value_cutoff=0.05,
                      adj.P.Val_cutoff=0.05, 
                      dPDUI_cutoff=0.3, 
                      PDUI_logFC_cutoff){
    if(missing(gp1) || missing(gp2)){
        testRes <- res$testRes
        testRes <- testRes[!(is.na(testRes[, "P.Value"]) | 
                                 is.na(testRes[, "adj.P.Val"])), ]
        if(!"dPDUI" %in% colnames(testRes)){
            message("dPDUI_cutoff will be ignored")
        }else{
            testRes <- testRes[abs(testRes[, "dPDUI"])>=dPDUI_cutoff, , drop=FALSE]
        }
        if(!missing(PDUI_logFC_cutoff)){
            if(!"logFC" %in% colnames(testRes)){
                message("PDUI_logFC_cutoff will be ignored")
            }else{
                testRes <- testRes[abs(testRes[, "logFC"])>=PDUI_logFC_cutoff, , drop=FALSE]
            }
        }
        testRes <- testRes[testRes[, "adj.P.Val"]<=adj.P.Val_cutoff &
                               testRes[, "P.Value"]<=P.Value_cutoff, , drop=FALSE]
        testRes <- testRes[!is.na(rownames(testRes)), ]
        if(nrow(testRes)==0) return(NULL)
        res <- res$usage[match(rownames(testRes), res$usage$transcript)]
        mcols(res) <- cbind(as.data.frame(mcols(res)), testRes)
        return(res)
    }
    ## for each group, at least in one group, more than half of the long form should greater than background_coverage_threshold
    ## for both group, at leat in one group, more than half of the short form should greater than background_coverage_threshold
    coln <- colnames(res$short)
    if(!all(c(gp1, gp2) %in% coln)){
        stop("gp1 and gp2 must be the tags of samples (colnames of res$short)")
    }
    short <- res$short
    long <- res$long
    short.abg <- short > background_coverage_threshold
    long.abg <- long > background_coverage_threshold
    short.abg.gp1 <- short.abg[, gp1, drop=FALSE]
    short.abg.gp2 <- short.abg[, gp2, drop=FALSE]
    long.abg.gp1 <- long.abg[, gp1, drop=FALSE]
    long.abg.gp2 <- long.abg[, gp2, drop=FALSE]
    halfLgp1 <- ceiling(length(gp1)/2)
    halfLgp2 <- ceiling(length(gp2)/2)
    short.abg.gp1 <- rowSums(short.abg.gp1)
    short.abg.gp2 <- rowSums(short.abg.gp2)
    long.abg.gp1 <- rowSums(long.abg.gp1)
    long.abg.gp2 <- rowSums(long.abg.gp2)
    long.abg <- long.abg.gp1 >= halfLgp1 | long.abg.gp2 >= halfLgp2
    short.abg <- short.abg.gp1 >= halfLgp1 | short.abg.gp2 >= halfLgp2
    abg.gp1 <- long.abg.gp1 >= halfLgp1 | short.abg.gp1 >= halfLgp1
    abg.gp2 <- long.abg.gp2 >= halfLgp2 | short.abg.gp2 >= halfLgp2
    abg <- abg.gp1 & abg.gp2 & long.abg & short.abg  
    
    PDUI.gp1 <- res$PDUI[, gp1, drop=FALSE]
    PDUI.gp2 <- res$PDUI[, gp2, drop=FALSE]
    PDUI.gp1 <- rowMeans(PDUI.gp1)
    PDUI.gp2 <- rowMeans(PDUI.gp2)
    if(!"dPDUI" %in% colnames(res$testRes)){
        message("dPDUI is calculated by gp2 - gp1.")
        dPDUI <- PDUI.gp2 - PDUI.gp1
    }else{
        dPDUI <- res$testRes[, "dPDUI"]
    }
    adPDUI <- abs(dPDUI) >= dPDUI_cutoff
    
    if(!missing(PDUI_logFC_cutoff)){
        if(!"logFC" %in% colnames(res$testRes)){
            logFC <- log2(PDUI.gp2 + .Machine$double.xmin) - 
                log2(PDUI.gp1 + .Machine$double.xmin)
        }else{
            logFC <- res$testRes[, "logFC"]
        }
        alFC <- abs(logFC) >= PDUI_logFC_cutoff
    }else{
        alFC <- rep(TRUE, length(abg))
    }
    
    res <- as(res, "GRanges")
    res$dPDUI <- dPDUI
    res$PASS <- res$adj.P.Val<=adj.P.Val_cutoff &
        res$P.Value<=P.Value_cutoff &
        abg & adPDUI & alFC
    res
}