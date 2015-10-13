testUsage <- function(CPsites, coverage, genome, utr3, BPPARAM=NULL, 
                      method=c("limma", 
                                "fisher.exact", 
                                "singleSample", 
                                "singleGroup"),
                      normalize=c("none", "quantiles", "quantiles.robust",
                               "mean", "median"),
                      design, contrast.matrix, coef=1, robust=FALSE, ..., 
                      gp1, gp2){
    method <- match.arg(method)
    normalize <- match.arg(normalize)
    res <- NULL
    eset <- list()
    if(method!="singleSample"){
        eset <- getUTR3eSet(CPsites, coverage, genome, utr3, normalize, BPPARAM=BPPARAM)
        if(method=="limma"){
            if(missing(design)||missing(contrast.matrix)){
                stop("desing and contrast.matrix is required.")
            }
            res <- limmaAnalyze(eset, design, contrast.matrix, coef=coef, robust=FALSE, ...)
        }
        if(method=="fisher.exact"){
            if(missing(gp1)||missing(gp2)){
                stop("gp1 and gp2 is required.")
            }
            res <- 
                fisher.exact.test(eset, gp1=gp1, gp2=gp2)
        }
        if(method=="singleGroup"){
            res <- singleGroupAnalyze(eset)
        }
    }else{
        eset <- getUTR3eSet(CPsites, coverage, genome, utr3, BPPARAM=BPPARAM, singleSample=TRUE)
        res <- singleSampleAnalyze(eset)
    }
    eset$testRes <- res
    return(eset)
}