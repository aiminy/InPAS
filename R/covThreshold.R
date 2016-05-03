intronRegion <- function(txdb){
    if(class(txdb)!="TxDb") stop("txdb must be an object of TxDb")
    exons <- ranges(GeneRegionTrack(txdb))
    exons <- reduce(exons)
    genes <- genes(txdb)
    genes <- reduce(genes)
    dis <- disjoin(c(genes, exons))
    ol <- findOverlaps(dis, exons)
    intronRegion <- dis[-unique(queryHits(ol))]
}

## auto determin the long_coverage_threshold by 
## non_zero intergenicRegion coverage
## auto determin the short_coverage_threshold by 
## quatile of non_zero intragenicRegion coverage
covThreshold <- function(coverage, genome, txdb, utr3, 
                         chr="chr1", hugeData, groupList){
    totalCov <- totalCoverage(coverage, genome, hugeData, groupList)
    chr1totCov <- lapply(totalCov, "[[", chr)
    N <- length(chr1totCov)
    if(N>1){
        for(i in 2:N){
            chr1totCov[[1]] <- chr1totCov[[1]] + chr1totCov[[i]]
        }
    }
    chr1totCov <- chr1totCov[[1]]
    intronRegion <- intronRegion(txdb)
    intronRegion <- intronRegion[seqnames(intronRegion)==chr]
    utr3Chr1Region <- utr3[seqnames(utr3)==chr]
    utr3Chr1Region <- utr3Chr1Region[utr3Chr1Region$id=="utr3"]
    covBg<-function(.cvg, start, end){
        view <- Views(.cvg, start, end)#
        view <- viewApply(view, function(.ele) as.integer(.ele))
        view <- unlist(view)
        view <- view[view>=N]
        floor(quantile(view)[2])
    }
    long_coverage_threshold <- covBg(chr1totCov, 
                                     start(intronRegion), 
                                     end(intronRegion))
    short_coverage_threshold <- covBg(chr1totCov, 
                                      start(utr3Chr1Region), 
                                      end(utr3Chr1Region))
    ceiling(c(long_coverage_threshold/N, short_coverage_threshold/N))
}
