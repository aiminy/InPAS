limmaAnalyze <- function(UTR3eset, design, contrast.matrix){
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
    fit <- lmFit(UTR3eset$PDUI.log2, design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit1)
    return(fit2)
}