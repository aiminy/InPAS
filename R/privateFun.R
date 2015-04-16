trimSeqnames <- function(genome){
    if(class(genome)!="BSgenome") stop("genome must be an object of BSgenome")
    seqnames <- seqnames(genome)
    seqnames <- seqnames[grepl("^chr[0-9XY]+$", seqnames)]
}

seqLen <- function(genome){
    if(class(genome)!="BSgenome") stop("genome must be an object of BSgenome")
    seqnames <- trimSeqnames(genome)
    seqlengths(genome)[seqnames]
}

depthWeight <- function(coverage, hugeData, gp1=NA, gp2=NA){
    if(hugeData){
        if((is.character(gp1)) && (is.character(gp2))){
            depth <- numeric(2)
            names(depth) <- c("gp1", "gp2")
            for(i in 1:length(coverage)){
                cvg <- NULL
                if(names(coverage)[i] %in% gp1){
                    j <- "gp1"
                }else{
                    if(names(coverage)[i] %in% gp2){
                        j <- "gp2"
                    }else{
                        j <- NA
                    }
                }
                if(!is.na(j)){
                    load(coverage[[i]])
                    d <- sapply(cvg, function(.cvg) {
                        sum(as.double(runValue(.cvg)) * runLength(.cvg))
                    })
                    depth[j] <- sum(c(d), depth[j])
                    rm(cvg)
                }
            }
        }else{
            depth <- numeric(length(coverage))
            for(i in 1:length(coverage)){
                cvg <- NULL
                load(coverage[[i]])
                d <- sapply(cvg, function(.cvg) {
                    sum(as.double(runValue(.cvg)) * runLength(.cvg))
                })
                depth[i] <- sum(d)
                rm(cvg)
            }
        }
    }else{
        depth <- sapply(coverage, function(cvg){
            d <- sapply(cvg, function(.cvg) {
                sum(as.double(runValue(.cvg)) * runLength(.cvg))
            })
            sum(d)
        })
    }
    depth.weight <- depth / mean(depth)
    if(hugeData && (is.character(gp1)) && (is.character(gp2))){
        names(depth.weight) <- c("gp1", "gp2")
    }else{
        names(depth.weight) <- names(coverage)
    }
    
    depth.weight
}

sortGR <- function(.ele){
    str <- as.character(strand(.ele))[1] == "+"
    if(str){ ## sort the element by position
        .ele <- .ele[order(start(.ele))]
    }else{
        .ele <- .ele[order(-end(.ele))]
    }
    .ele
}

getCov <- function(bedgraph, genome, seqLen){
    seqnames.bedfile <- 
        read.delim(bedgraph, header=FALSE, comment.char="#", 
                   colClasses=c(NA, "NULL", "NULL", "NULL"))[,1]
    seqnames <- sort(intersect(levels(seqnames.bedfile), trimSeqnames(genome)))
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
                         chr="chr1", hugeData, gp1, gp2){
    totalCov <- totalCoverage(coverage, genome, hugeData, gp1, gp2)
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


## calcultate p value:
## T=abs(mean(sqrt(short_coverage1))/mean(sqrt(long_coverage2)) 
##    - mean(sqrt(long_coverage1))/mean(sqrt(long_coverage2))) 
## sqrt also could be log2, sqrt --> poisson distribution
## log2 --> power law
## p.value = n(Ti<T)/n.total(T)

valley <- function(x, ss, se, n=1, savedID=NA, filterByPval=TRUE){
    pos <- which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == 2) + ss - 1
    if(length(pos)>0 && filterByPval){
        pos.pos <- 
            which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == -2) + ss - 1
        pos.pos <- sort(c(pos, pos.pos, ss, se))
        w <- which(pos.pos %in% pos)
        space=10
        y <- mapply(function(a, b, p) {
            c <- unique(sort(seq(a, b)))
            c <- c[c!=p & c>ss & c<se]
            d <- x[c]
        }, pos.pos[w-1] - space, pos.pos[w-1] + space, pos, SIMPLIFY=FALSE)
        mu <- sapply(y, mean, na.rm=TRUE)
        sigma <- sapply(y, sd, na.rm=TRUE)
        z <- (x[pos]-mu)/sigma
        p <- 2*pnorm(-abs(z))
        pos <- pos[p < 0.001]
    }
    #pos <- unique(c(ss, pos, se))
    pos <- pos[!is.na(x[pos])]
    tobeadd <-ifelse(x[ss] < x[se], ss, se)
    if(length(pos)>0){
        w <- x[tobeadd] < min(x[pos])
        if(length(w)>0 && !is.na(w) && is.logical(w)){
            if(w) pos <- unique(c(pos, tobeadd))
        }
    }else{
        pos <- tobeadd
    }
    if(!is.na(savedID) && savedID>ss && savedID<se) 
        pos <- unique(c(pos, savedID))
    pos <- pos[x[pos]<=quantile(x[ss:se], probs=.5)]
    if(n==-1 || n>length(pos)) n <- length(pos)
    if(length(pos)<1){
        return(pos)
    }else{
        pos <- pos[order(x[pos])][1:n]
    }
}