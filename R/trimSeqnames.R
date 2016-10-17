trimSeqnames <- function(genome){
    if(class(genome)!="BSgenome") stop("genome must be an object of BSgenome")
    seqnames <- seqnames(genome)
    seqnames <- seqnames[grepl("^chr[0-9XY]+$", seqnames)]
}