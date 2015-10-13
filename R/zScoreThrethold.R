zScoreThrethold <- function(background, introns, totalCov, utr3, z=2){
    if(background=="same_as_long_coverage_threshold")
        return(0)
    if(length(introns)<1)
        return(0)
    
    curr_range <- utr3[utr3$id=="utr3"]
    
    seqnames <- unlist(sapply(totalCov, names))
    seqnames <- table(seqnames)
    seqnames <- names(seqnames)[seqnames==length(totalCov)]
    seqnames <- sort(intersect(as.character(seqnames(curr_range)), 
                               seqnames))
    
    background <- switch(background, 
                         "1K"=1000,
                         "5K"=5000,
                         "10K"=10000,
                         "50K"=50000,
                         1000)
    
    curr_str <- as.character(strand(curr_range)) == "+"
    start(curr_range) <- ifelse(curr_str, 
                                start(curr_range)-background-1,
                                end(curr_range))
    width(curr_range) <- background
    
    curr_range <- curr_range[seqnames(curr_range) %in% seqnames]
    introns <- introns[seqnames(introns) %in% seqnames]
    ol <- findOverlaps(curr_range, introns)
    if(length(ol)<1)
        return(0)
    introns <- introns[subjectHits(ol)]
    introns$id <- names(curr_range)[queryHits(ol)]
    
    cov <- list()
    for(seq in seqnames){
        cov[[seq]] <- list()
        for(n in names(totalCov)){
            cov[[seq]][[n]] <- totalCov[[n]][[seq]]
        }
    }
    
    introns.s <- split(introns, as.character(seqnames(introns)))
    cvg <- mapply(function(.intron, .cov){
        .cvg <- 
            lapply(.cov, function(.cv){
                .cv <- Views(.cv, 
                             start=start(.intron), 
                             end=end(.intron))
                .cv <- viewApply(.cv, as.integer, simplify=FALSE)
            })
        .cvg1 <- vector("list", length(.cvg))
        for(i in 1:length(.cvg[[1]])){
            .cvg1[[i]] <- list()
        }
        for(i in names(.cvg)){
            for(j in 1:length(.cvg[[1]])){
                .cvg1[[j]][[i]] <- .cvg[[i]][[j]]
            }
        }
        .cvg1 <- lapply(.cvg1, function(.ele){
            colSums(do.call(rbind, .ele))
        })
        .intron$cvg <- .cvg1
        .intron
    }, introns.s[seqnames], cov[seqnames], SIMPLIFY=FALSE)
    
    cvg <- unlist(GRangesList(cvg))
    cvg <- split(cvg$cvg, cvg$id)
    cvg <- sapply(cvg, function(.ele){
        .ele <- unlist(.ele)
        mu <- mean(.ele, na.rm=TRUE)
        std <- sd(.ele, na.rm=TRUE)
        ##Z-score = (x-mu)/std
        ##pval = 2*pnorm(-abs(z))
        ##Z-score = qnorm(1 - (pval/2)) 
        ##Z-score should greate than 2 for p=.05
        ## z < (x-mu)/std
        ## z*std+mu < x
        z*std + mu
    })
    cvg <- cvg[match(names(curr_range), names(cvg))]
    names(cvg) <- names(curr_range)
    cvg
}