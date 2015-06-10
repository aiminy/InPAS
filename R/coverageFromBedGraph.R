coverageFromBedGraph <- function(bedgraphs, tags, 
                                 genome, hugeData=FALSE, ...){
    if(missing(genome))
        stop("genome is required.")
    if(class(genome)!="BSgenome")
        stop("genome must be an object of BSgenome.")
    if(missing(tags) || missing(bedgraphs))
        stop("tags and bedgraphs are required.")
    if(length(tags)!=length(bedgraphs)){
        stop("length of tags and bedgraphs should be identical")
    }
    if(class(tags)!="character")
        stop("tags must be a character vector")
    if(any(duplicated(tags)))
        stop("There are duplicated tags")
    if(any(!file.exists(bedgraphs)))
        stop("Not all bedgraphs exist")
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