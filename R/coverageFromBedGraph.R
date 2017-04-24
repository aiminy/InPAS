coverageFromBedGraph <- function(bedgraphs, tags, 
                                 genome, hugeData=FALSE, 
                                 BPPARAM=NULL, ...){
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
    null <- sapply(bedgraphs, function(.ele){
        ## check the bedgraphs, can not be empty.
        read.delim(.ele, header=FALSE, 
                   comment.char="#", nrows=5,
                   colClasses=c("factor", "NULL", "NULL", "NULL"))
    })
    seqLen <- seqLen(genome)
    ## get coverage for all inputs
    names(bedgraphs) <- tags
    x <- 1:length(bedgraphs)
    y <- split(x, ceiling(x/10))
    if(hugeData){
        coverage <- bedgraphs
        for(i in 1:length(y)){
            if(!is.null(BPPARAM)){
                cv <- bptry(bplapply(bedgraphs[y[[i]]], function(.ele){
                    cvg <- getCov(.ele, genome, seqLen)
                    filename <- tempfile(...)
                    save(list="cvg", file=filename)
                    filename
                }, BPPARAM=BPPARAM))
                while(!all(bpok(cv))){
                    cv <- bptry(bplapply(bedgraphs[y[[i]]], function(.ele){
                        cvg <- getCov(.ele, genome, seqLen)
                        filename <- tempfile(...)
                        save(list="cvg", file=filename)
                        filename
                    }, BPREDO=cv, BPPARAM=BPPARAM))
                }
            }else{
                cv <- lapply(bedgraphs[y[[i]]], function(.ele){
                    cvg <- getCov(.ele, genome, seqLen)
                    filename <- tempfile(...)
                    save(list="cvg", file=filename)
                    filename
                })
            }
            coverage[names(cv)] <- unlist(cv)
        }
        coverage <- as.list(coverage[tags])
    }else{
        coverage <- list()
        for(i in 1:length(y)){
            if(!is.null(BPPARAM)){
                cv <- bplapply(bedgraphs, getCov, 
                                     genome=genome, seqLen=seqLen,
                                     BPPARAM=BPPARAM)
            }else{
                cv <- lapply(bedgraphs, getCov, genome=genome, seqLen=seqLen)
            }
            coverage <- c(coverage, cv)
        }
        coverage <- coverage[tags]
    }
    coverage
}