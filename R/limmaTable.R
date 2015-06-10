limmaTable <- function(fit, UTR3eset, coef=1, ...){
    if(class(UTR3eset)!="UTR3eSet"){
        stop("UTR3eset must be an object of UTR3eSet")
    }
    res <- topTable(fit, coef=coef, number=nrow(fit), ...)
    res <- res[UTR3eset$usage$transcript, , drop=FALSE]
    res
}