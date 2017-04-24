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
    
    ex <- diffSplice(fit, geneid=txid, exonid=formid, 
                     robust=robust, verbose=FALSE)
    ts <- topSplice(ex, coef=coef, test="simes", number=nrow(ex))
    p <- ts[match(rownames(short), ts$GeneID), "P.Value"]
    BH <- ts[match(rownames(short), ts$GeneID), "FDR"]
    out <- cbind(P.Value=p, adj.P.Val=BH)
    rownames(out) <- rownames(short)
    out
}