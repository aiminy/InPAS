limmaAnalyze <- function(UTR3eset, design, contrast.matrix, coef=1, robust=FALSE, ...){
    if(missing(design)||missing(contrast.matrix)){
        stop("desing and contrast.matrix is required.")
    }
    if(class(design)!="matrix"){
        stop("design must be an design matrix")
    }
    if(class(contrast.matrix)!="matrix"){
        stop("contrast.matrix must be an object of matrix")
    }
    if(class(UTR3eset)!="UTR3eSet"){
        stop("UTR3eset must be an object of UTR3eSet")
    }
    short <- UTR3eset@short
    long <- UTR3eset@long
    if(!identical(rownames(short), rownames(long))){
        stop("the rownames are not identical for short form and long form")
    }
    txid <- rep(rownames(short), 2)
    formid <- c(paste(rownames(short), "short", sep=":"), 
                paste(rownames(long), "long", sep=":"))
    y <- rbind(short, long)
    rownames(y) <- formid
    lib.size <- colSums(y)
    y <- log2(y+0.001)
    y <- normalizeBetweenArrays(y)
    fit <- lmFit(y, design, ...)
    fit <- contrasts.fit(fit, contrast.matrix)
    
    genecolname <- "txid"
    exoncolname <- "formid"
    form.tx <- data.frame(formid=formid, txid=txid)
    
    ## sort by txid
    o <- order(txid, formid)
    txid <- txid[o]
    form.tx <- form.tx[o, , drop=FALSE]
    
    form.coefficients <- fit$coefficients[o,,drop=FALSE]
    form.stdev.unscaled <- fit$stdev.unscaled[o,,drop=FALSE]
    form.df.residual <- fit$df.residual[o]
    form.s2 <- fit$sigma[o]^2
    
    ##Count forms and get transcript-wise variances
    form.stat <- cbind(1,form.df.residual,form.s2)
    tx.sum <- rowsum(form.stat,txid,reorder=FALSE)
    tx.nforms <- tx.sum[,1]
    if(!all(tx.nforms==2)){
        stop("Bug when check the length of forms")
    }
    tx.df.residual <- tx.sum[,2]
    tx.s2 <- tx.sum[,3] / tx.sum[,1]
    squeeze <- squeezeVar(var=tx.s2, df=tx.df.residual, robust=robust)
    
    tx.df.test <- tx.nforms-1
    tx.df.total <- tx.df.residual+squeeze$df.prior
    tx.df.total <- pmin(tx.df.total,sum(tx.df.residual))
    tx.s2.post <- squeeze$var.post
    
    #     transcript-wise betas
    u2 <- 1/form.stdev.unscaled^2
    u2.rowsum <- rowsum(u2, txid,reorder=FALSE)
    tx.betabar <- rowsum(form.coefficients*u2,txid,reorder=FALSE) / u2.rowsum
    
    #    T-statistics for form-level tests
    ntx <- nrow(short)
    g <- rep(1:ntx,times=tx.nforms)
    form.coefficients <- form.coefficients-tx.betabar[g,,drop=FALSE]
    form.t <- form.coefficients / form.stdev.unscaled / sqrt(tx.s2.post[g])
    tx.F <- rowsum(form.t^2,txid,reorder=FALSE) / tx.df.test
    form.1mleverage <- 1 - (u2 / u2.rowsum[g,,drop=FALSE])
    form.coefficients <- form.coefficients / form.1mleverage
    form.t <- form.t / sqrt(form.1mleverage)
    form.p.value <- 2 * pt(abs(form.t), df=tx.df.total[g], lower.tail=FALSE)
    tx.F.p.value <- pf(tx.F, df1=tx.df.test, df2=tx.df.total, lower.tail=FALSE)
    
    p <- tx.F.p.value[, coef]
    BH <- p.adjust(p, method="BH")
    out <- cbind(P.Value=p, adj.P.Val=BH)
}