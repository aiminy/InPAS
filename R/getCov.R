getCov <- function(bedgraph, genome, seqLen){
    seqnames.bedfile <- 
        read.delim(bedgraph, header=FALSE, comment.char="#", 
                   colClasses=c("factor", "NULL", "NULL", "NULL"))[,1]
    seqnames <- trimSeqnames(genome)
    seqStyle <- seqlevelsStyle(genome) ## should be UCSC
    seqStyle.bed <- seqlevelsStyle(levels(seqnames.bedfile))
    if(!any(seqStyle.bed==seqStyle)){
        message("seqlevelsStyle of genome is different from bedgraph file.")
        ## convert to seqStyle
        levels <- 
            mapSeqlevels(levels(seqnames.bedfile), seqStyle)
        if(class(levels)!="character"){
            id <- apply(levels, 1, function(.ele) 
                sum(seqnames %in% .ele))
            levels <- levels[id==max(id),][1, ]
        }
        levels[is.na(levels)] <- names(levels)[is.na(levels)]
        levels(seqnames.bedfile) <- levels
    }
    seqnames <- sort(intersect(levels(seqnames.bedfile), seqnames))
    if(length(seqnames)<1){
        stop(paste("there is no intersect seqname in", bedgraph, "and genome"))
    }
    summaryFunction <- function(seqname){
        seqL <- seqLen[seqname]
        lines2read <- Rle(c(FALSE, seqnames.bedfile==seqname))
        true <- which(runValue(lines2read))
        skip <- runLength(lines2read)[true-1]
        skip[1] <- skip[1] - 1
        nrow <- runLength(lines2read)[true]
        
        con <- file(bedgraph)
        open(con)
        dat <- read.table(con, nrows=nrow[1], header=FALSE, sep="\t", 
                          skip=skip[1], 
                          colClasses=c("NULL", "integer", "integer", NA))
        if(length(true)>1){
            for(i in 2:length(true)){
                lines <- 
                    read.table(con, nrows=nrow[i], header=FALSE, sep="\t", 
                               skip=skip[i], 
                               colClasses=c("NULL", "integer", "integer", NA))
                dat <- rbind(dat, lines)
            }
        }
        close(con)
        
        dat[,1] <- dat[,1]+1
        gaps <- 
            as.data.frame(gaps(IRanges(dat[,1], dat[,2]), start=1, end=seqL))
        if(nrow(gaps)>0){
            gaps <- cbind(gaps[,1:2], V4=0)
            colnames(gaps) <- colnames(dat)
            dat <- rbind(dat, gaps)
        }
        dat <- dat[order(dat[,1]),]
        Rle(dat[,3], dat[,2]-dat[,1]+1)
    }
    cov <- lapply(seqnames, summaryFunction)
    names(cov) <- seqnames
    cov
}