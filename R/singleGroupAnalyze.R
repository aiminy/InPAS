singleGroupAnalyze <- function(UTR3eset){
    data.long <- UTR3eset$long
    data.short <- UTR3eset$short
    data <- log2(cbind(data.long, data.short)+.Machine$double.xmin)
    treatments <- cbind(long=c(rep(c(1,0),c(ncol(data.long), 
                                            ncol(data.short)))), 
                        short=c(rep(c(0,1), c(ncol(data.long), 
                                              ncol(data.short)))))
    design <- model.matrix(~-1+treatments)
    colnames(design) <- c("long", "short")
    fit <- lmFit(data, design)
    contrast.matrix<-makeContrasts(contrasts="long-short",
                                   levels=design)
    fit <- contrasts.fit(fit,contrast.matrix)
    fit <- eBayes(fit)
    data <- topTable(fit, number=nrow(fit), sort.by="none")
    P.Value <- data$P.Value
    adj.P.Val <- p.adjust(P.Value, method="BH")
    long <- rowMeans(data.long)
    short <- rowMeans(data.short)
    PDUI <- long/c(long+short)
    cbind(short.mean=short, 
          long.mean=long,
          PDUI, P.Value, adj.P.Val)
}