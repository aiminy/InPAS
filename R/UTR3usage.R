UTR3usage <- function(CPsites, coverage, hugeData, BPPARAM=NULL, phmm=FALSE){
    utr3.selected.shorten.UTR <- 
        split(CPsites, as.character(seqnames(CPsites)))
    seqnames <- names(utr3.selected.shorten.UTR)
    
    ##covert coverage from list of Rle to list of matrix
    ##step1 prepare the short and long utr3 regions
    utr3.trans.shorten.UTR <- split(CPsites, CPsites$transcript)
    if(!is.null(BPPARAM)){
        utr3.regions <- bplapply(utr3.trans.shorten.UTR, 
                                 getUTR3region, BPPARAM=BPPARAM)
    }else{
        utr3.regions <- lapply(utr3.trans.shorten.UTR, getUTR3region)
    }
    utr3.regions <- unlist(GRangesList(utr3.regions))
    ##step2, calculate coverage and merge into a matrix for long and short
    utr3.regions.chr <- split(utr3.regions, as.character(seqnames(utr3.regions)))
    utr3.regions.chr <- utr3.regions.chr[seqnames]
    if(!is.null(BPPARAM)){
        utr3.regions.chr <- bplapply(seqnames, get.regions.coverage, 
                                     BPPARAM=BPPARAM, 
                                     utr3.regions.chr=utr3.regions.chr, 
                                     hugeData=hugeData, 
                                     coverage=coverage, 
                                     phmm=phmm)
    }else{
        utr3.regions.chr <- lapply(seqnames, get.regions.coverage, 
                                   utr3.regions.chr=utr3.regions.chr, 
                                   hugeData=hugeData, 
                                   coverage=coverage, 
                                   phmm=phmm)
    }
    cov.ranges <- unlist(GRangesList(utr3.regions.chr))
}