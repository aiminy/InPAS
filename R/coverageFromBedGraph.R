coverageFromBedGraph <- function(bedgraphs, tags, 
                                 genome, hugeData=FALSE, ...){
    seqLen <- seqLen(genome)
    ## get coverage for all inputs
    if(hugeData){
        coverage <- list()
        for(i in 1:length(bedgraphs)){
            cvg <- getCov(bedgraphs[i], genome, seqLen)
            filename <- tempfile(...)
            save(list="cvg", file=filename)
            coverage[[tags[i]]] <- filename
        }
    }else{
        coverage <- lapply(bedgraphs, getCov, genome=genome, seqLen=seqLen)
        names(coverage) <- tags
    }
    coverage
}