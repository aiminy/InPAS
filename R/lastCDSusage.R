lastCDSusage <- function(CDS, coverage, hugeData, BPPARAM=NULL, phmm=FALSE){
    CDS.regions.chr <- split(CDS, as.character(seqnames(CDS)))
    seqnames <- names(CDS.regions.chr)
    if(!is.null(BPPARAM)){
        CDS.regions.chr <- bplapply(seqnames, get.regions.coverage, 
                                    BPPARAM=BPPARAM,
                                    utr3.regions.chr=CDS.regions.chr,
                                    hugeData=hugeData,
                                    coverage=coverage,
                                    phmm=phmm)
    }else{
        CDS.regions.chr <- lapply(seqnames, get.regions.coverage, 
                                  utr3.regions.chr=CDS.regions.chr,
                                  hugeData=hugeData,
                                  coverage=coverage,
                                  phmm=phmm)
    }
    cov.ranges <- unlist(GRangesList(CDS.regions.chr))
}