fisher.exact.test <- function(UTR3eset, gp1, gp2){
    short <- UTR3eset$short
    long <- UTR3eset$long
    if(!all(c(gp1,gp2) %in% colnames(UTR3eset$short))) 
        stop("gp1 and gp2 must be same as your samples")
    short.gp1 <- short[, gp1, drop=FALSE]
    short.gp2 <- short[, gp2, drop=FALSE]
    long.gp1 <- long[, gp1, drop=FALSE]
    long.gp2 <- long[, gp2, drop=FALSE]
    short.mean.gp1 <- rowMeans(short.gp1)
    short.mean.gp2 <- rowMeans(short.gp2)
    long.mean.gp1 <- rowMeans(long.gp1)
    long.mean.gp2 <- rowMeans(long.gp2)
    mat <- cbind(short.mean.gp1,
                 short.mean.gp2,
                 long.mean.gp1,
                 long.mean.gp2)
    mat[is.na(mat)] <- 0
    P.Value <- apply(mat, 1, function(.d) 
        fisher.test(matrix(floor(.d), ncol=2))$p.value)
    adj.P.Val <- p.adjust(P.Value, method="BH")
    PDUI.gp1 <- long.mean.gp1 / (long.mean.gp1 + short.mean.gp1)
    PDUI.gp2 <- long.mean.gp2 / (long.mean.gp2 + short.mean.gp2)
    dPDUI <- PDUI.gp1 - PDUI.gp2
    logFC <- log2(PDUI.gp1+.Machine$double.xmin) - 
        log2(PDUI.gp2+.Machine$double.xmin)
    cbind(short.mean.gp1, long.mean.gp1, 
          short.mean.gp2, long.mean.gp2,
          PDUI.gp1, PDUI.gp2, dPDUI, logFC, 
          P.Value, adj.P.Val)
}