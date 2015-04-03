totalCoverage <- function(coverage, genome, hugeData, gp1, gp2){
    if(!hugeData) return(coverage)
    ## calculate total coverage
    cov <- list()
    if(is.character(gp1) && is.character(gp2)){
        cov[["gp1"]] <- list()
        cov[["gp2"]] <- list()
        seqnames <- trimSeqnames(genome)
        seqLen <- seqLen(genome)
        for(i in 1:2){
            for(s in seqnames){
                cov[[i]][[s]] <- Rle(0, seqLen[s])
            }
        }
        
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            gp <- names(coverage)[i]
            gp <- ifelse(gp %in% gp1, "gp1", ifelse(gp %in% gp2, "gp2", NA))
            if(!is.na(gp)){
                for(s in seqnames){
                    if(s %in% names(cvg)){
                        cov[[gp]][[s]] <- cov[[gp]][[s]] + cvg[[s]]
                    }
                }
            }
            rm(cvg)
        }
        idx <- rep(FALSE, length(seqnames))
        names(idx) <- seqnames
        for(i in 1:2){
            for(s in seqnames){
                if(nrun(cov[[i]][[s]])!=1 || runValue(cov[[i]][[s]][1])!=0) 
                    idx[s] <- TRUE
            }
        }
        for(i in 1:2){
            cov[[i]] <- cov[[i]][idx]
        }
    }else{
        seqnames <- trimSeqnames(genome)
        seqLen <- seqLen(genome)
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            gp <- names(coverage)[i]
            cov[[gp]] <- list()
            for(s in seqnames){
                if(s %in% names(cvg)){
                    cov[[gp]][[s]] <- cvg[[s]]
                }
            }
            rm(cvg)
        }
    }
    cov
}


compensation <- function(view, comp, start, end){
    mapply(function(.ele, .s, .e){
        .ele * comp[.s:.e]
    }, view, start, end, SIMPLIFY=FALSE)
}

UTR3TotalCoverage <- function(utr3, totalCov, 
                              gcCompensation=NA, 
                              mappabilityCompensation=NA, 
                              FFT=FALSE, 
                              fft.sm.power=20){
    utr3.s <- split(utr3, as.character(seqnames(utr3)))
    ##smooth methods
    fft.smooth <- function(sn, p){#
        sn.fft <- fft(sn)#
        sn.fft[p:length(sn.fft)] <- 0+0i#
        sn.ifft = fft(sn.fft, inverse = TRUE)/length(sn.fft)#
        Re(sn.ifft)#
    }
    seqnames <- unlist(sapply(totalCov, names))
    seqnames <- table(seqnames)
    seqnames <- names(seqnames)[seqnames==length(totalCov)]
    seqnames <- sort(intersect(names(utr3.s), seqnames))
    cov <- list()
    for(seq in seqnames){
        cov[[seq]] <- list()
        for(n in names(totalCov)){
            cov[[seq]][[n]] <- totalCov[[n]][[seq]]
        }
    }
    mapply(function(.utr, .cvg, gcComp, mapComp){ ##memory consume
        strand <- as.character(strand(.utr))
        start <- start(.utr)
        end <- end(.utr)
        maxEnd <- max(end)
        for(i in 1:length(.cvg)){
            if(maxEnd>length(.cvg[[i]])){
                .cvg[[i]] <- 
                    append(.cvg[[i]], rep(0, maxEnd - length(.cvg[[i]])+1))
            }
        }
        
        view <- lapply(.cvg, Views, start=start, end=end, names=names(.utr))
        view <- sapply(view, function(.view){
            viewApply(.view, as.integer)
        })
        if(class(view[1,])=="list"){
            view <- apply(view, 1, function(.ele) do.call(cbind, .ele))
        }else{
            view <- list(view)
            names(view) <- names(.utr)
        }
        if(length(view)>0) view <- view[sapply(view, sum)>0]
        if(!is.na(gcComp[1]) && length(gcComp)==length(.cvg[[1]])){
            view <- compensation(view, gcComp, start, end)
        }
        if(!is.na(mapComp[1]) && length(mapComp)==length(.cvg[[1]])){
            view <- compensation(view, mapComp, start, end)
        }
        if(FFT){
            lapply(view, function(.ele) 
                apply(.ele, 2, fft.smooth, p=fft.sm.power))
        }else{
            view
        }
    }, utr3.s[seqnames], cov[seqnames], 
    gcCompensation[seqnames], 
    mappabilityCompensation[seqnames], 
    SIMPLIFY=FALSE)
}

CPsites <- function(coverage, gp1, gp2=NULL, genome, utr3, window_size=100, 
                    search_point_START=50, search_point_END=NA, 
                    cutStart=100, cutEnd=0, search_distal_polyA_end=FALSE, 
                    coverage_threshold=5, long_coverage_threshold=2, 
                    gcCompensation=NA, mappabilityCompensation=NA, 
                    FFT=FALSE, fft.sm.power=20, 
                    PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, 
                    shift_range=25, BPPARAM=NULL){
    if(!all(c(gp1,gp2) %in% names(coverage))) 
        stop("gp1 and gp2 must be in names of coverage")
    if(!is.na(PolyA_PWM)[1]){
        if(class(PolyA_PWM)!="matrix") stop("PolyA_PWM must be matrix")
        if(any(rownames(PolyA_PWM)!=c("A", "C", "G", "T"))) 
            stop("rownames of PolyA_PWM must be c(\"A\", \"C\", \"G\", \"T\")")
    }
    if(missing(coverage) || missing(genome) || missing(utr3))
        stop("coverage, genome and utr3 is required.")
    if(class(genome)!="BSgenome")
        stop("genome must be an object of BSgenome.")
    if(class(utr3)!="GRanges" | 
           !all(utr3$id %in% c("utr3", "next.exon.gap", "CDS"))){
        stop("utr3 must be output of function of utr3Annotation")
    }
    utr3 <- utr3[utr3$id!="CDS"]
    MINSIZE <- 10
    hugeData <- class(coverage[[1]])=="character"
    depth.weight <- depthWeight(coverage, hugeData, gp1, gp2)
    totalCov <- totalCoverage(coverage, genome, hugeData, gp1, gp2)
    utr3TotalCov <- UTR3TotalCoverage(utr3, totalCov, 
                                      gcCompensation, mappabilityCompensation, 
                                      FFT=FFT, fft.sm.power=fft.sm.power)
    if(!is.null(BPPARAM)){
        shorten_UTR_estimation <- 
            bplapply(utr3TotalCov, CPsite_estimation,
                     BPPARAM=BPPARAM, utr3=utr3,
                     MINSIZE=MINSIZE, 
                     window_size=window_size, 
                     search_point_START=search_point_START,
                     search_point_END=search_point_END, 
                     cutStart=cutStart, cutEnd=cutEnd, 
                     search_distal_polyA_end=search_distal_polyA_end, 
                     coverage_threshold=coverage_threshold, 
                     long_coverage_threshold=long_coverage_threshold, 
                     PolyA_PWM=PolyA_PWM, 
                     classifier=classifier, 
                     classifier_cutoff=classifier_cutoff, 
                     shift_range=shift_range, 
                     depth.weight=depth.weight, 
                     genome=genome)
    }else{
        shorten_UTR_estimation <- 
            lapply(utr3TotalCov, CPsite_estimation, 
                   utr3=utr3, MINSIZE=MINSIZE, 
                   window_size=window_size, 
                   search_point_START=search_point_START, 
                   search_point_END=search_point_END, 
                   cutStart=cutStart, cutEnd=cutEnd, 
                   search_distal_polyA_end=search_distal_polyA_end, 
                   coverage_threshold=coverage_threshold, 
                   long_coverage_threshold=long_coverage_threshold, 
                   PolyA_PWM=PolyA_PWM, 
                   classifier=classifier, classifier_cutoff=classifier_cutoff, 
                   shift_range=shift_range, 
                   depth.weight=depth.weight, 
                   genome=genome)
    }
    shorten_UTR_estimation <- 
        do.call(rbind, 
                shorten_UTR_estimation[!sapply(shorten_UTR_estimation, 
                                               is.null)])
    utr3.shorten.UTR <- utr3[utr3$id=="utr3"]
    utr3.shorten.UTR$id <- NULL
    utr3.shorten.UTR <- 
        utr3.shorten.UTR[utr3.shorten.UTR$transcript 
                         %in% rownames(shorten_UTR_estimation)]
    shorten_UTR_estimation <- 
        shorten_UTR_estimation[utr3.shorten.UTR$transcript, , drop=FALSE]
    utr3.shorten.UTR$fit_value <- unlist(shorten_UTR_estimation[,"fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <- 
        unlist(shorten_UTR_estimation[,"Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <- 
        unlist(shorten_UTR_estimation[,"Predicted_Distal_APA"])
    utr3.shorten.UTR$type <- unlist(shorten_UTR_estimation[,"type"])
    utr3.shorten.UTR$utr3start <- unlist(shorten_UTR_estimation[,"utr3start"])
    utr3.shorten.UTR$utr3end <- unlist(shorten_UTR_estimation[,"utr3end"])
    utr3.shorten.UTR <- 
        utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
}

CPsite_estimation <- function(chr.cov, utr3, MINSIZE, window_size, 
                              search_point_START, search_point_END, 
                              cutStart, cutEnd, search_distal_polyA_end, 
                              coverage_threshold, long_coverage_threshold, 
                              PolyA_PWM, classifier, classifier_cutoff, 
                              shift_range, depth.weight, genome){
    chr.cov <- chr.cov[sapply(chr.cov, sum)>0]
    if(length(chr.cov)==0){
        return(NULL)
    }
    curr_UTR.gr <- utr3[names(chr.cov)]
    utr3.utr <- curr_UTR.gr[curr_UTR.gr$id=="utr3"]
    utr3.gap <- curr_UTR.gr[curr_UTR.gr$id=="next.exon.gap"]
    co <- countOverlaps(utr3.gap, utr3.utr, maxgap=1, ignore.strand=TRUE)
    utr3.gap <- utr3.gap[co>1]
    curr_UTR.gr$conn_next_utr3 <- 
        curr_UTR.gr$transcript %in% utr3.gap$transcript
    curr_UTR <- split(curr_UTR.gr, curr_UTR.gr$transcript)
    conn_next_utr3 <- sapply(curr_UTR, function(.UTR){
        .UTR$conn_next_utr3[1]
    })
    chr.cov.merge <- lapply(curr_UTR, function(.UTR){
        .UTR <- .UTR[order(start(.UTR))]
        chr.utr3TotalCov <- chr.cov[names(.UTR)]
        chr.utr3TotalCov <- mapply(function(.covList, .start, .end, .property){
            #set names for each position 
            .posList <- .start:.end
            rownames(.covList) <- paste(.property, .posList, sep="_SEP_")
            .covList
        }, chr.utr3TotalCov, start(.UTR), end(.UTR), .UTR$id, SIMPLIFY=FALSE)
        chr.utr3TotalCov <- do.call(rbind, chr.utr3TotalCov)
        if(as.character(strand(.UTR))[1] == "-"){ ##reverse the negative strand
            chr.utr3TotalCov <- 
                chr.utr3TotalCov[rev(rownames(chr.utr3TotalCov)), , drop=FALSE]
        }
        if(!is.na(cutStart)){
            if(cutStart<1) cutStart <- floor(length(chr.utr3TotalCov)*cutStart)
            if(cutStart>0) 
                chr.utr3TotalCov <- 
                chr.utr3TotalCov[-(1:cutStart), , drop=FALSE]
        }
        
        chr.utr3TotalCov
    })
    ## chr.cov.merge should be a list with named numeric
    ## filter the coverage
    coverage_quality <- sapply(chr.cov.merge, function(.ele) {
        if(nrow(.ele[grepl("utr3_SEP_", rownames(.ele)), ,drop=FALSE])>
               MINSIZE){
            any(colMeans(.ele[1:min(nrow(.ele), 100), , drop=FALSE]) > 
                    coverage_threshold)
        }else{
            FALSE   
        }
    })
    
    chr.cov.merge <- chr.cov.merge[coverage_quality]
    conn_next_utr3 <- conn_next_utr3[coverage_quality]
    if(length(chr.cov.merge)>0){
        chr.abun <- mapply(function(chr.cov.merge.ele, 
                                    conn_next_utr, 
                                    curr_UTR.ele){
            chr.cov.merge.ele <- 
                t(t(chr.cov.merge.ele)/depth.weight[
                    colnames(chr.cov.merge.ele)])
            .ele <- rowSums(chr.cov.merge.ele)
            # if there are reads covered to the next.exon.gap
            # the proximal CP site should be the known-utr3-end, 
            # and distal site should be the end of gap
            next.exon.gap <- .ele[grepl("next.exon.gap", names(.ele))]
            # remove the gaps with 0 width > window_size
            next.exon.gap.rle <- Rle(next.exon.gap)
            id <- which(runLength(next.exon.gap.rle)>window_size & 
                            runValue(next.exon.gap.rle)==0)
            if(length(id)>=1){
                conn_next_utr <- FALSE
                id <- id[1] - 1
                if(id>0){
                    id <- sum(runLength(next.exon.gap.rle)[1:id])
                    next.exon.gap <- next.exon.gap[1:id]
                }else{
                    next.exon.gap <- numeric(0)   
                }
            }
            ## stop and any connected two windows different from 10 times, 
            ## to separate connected two UTRs or strange peaks.
            next.exon.gap.ids <- 
                rep(1:ceiling(length(next.exon.gap)/window_size), 
                    each=window_size)[1:length(next.exon.gap)]
            next.exon.gap.split <- split(next.exon.gap, next.exon.gap.ids)
            next.exon.gap.split <- sapply(next.exon.gap.split, mean)
            next.exon.gap.split.diff <- 
                log2(next.exon.gap.split[-1]+0.01) - 
                log2(next.exon.gap.split[-length(next.exon.gap.split)]+0.01)
            id <- which(abs(next.exon.gap.split.diff)>log2(10))
            if(length(id)>1){
                conn_next_utr <- FALSE
                id <- id[1] - 1
                if(id>0){
                    next.exon.gap <- next.exon.gap[1:(window_size*id)]
                }else{
                    next.exon.gap <- numeric(0)
                }
            }
            ## remove utr3---___---utr3, need to improve.
            if(conn_next_utr && length(next.exon.gap)>50){
                next.exon.gap <- removeUTR3__UTR3(next.exon.gap)
            }
            
            annotated.utr3 <- .ele[grepl("utr3", names(.ele))]
            utr3start <- as.numeric(gsub("^.*_SEP_", "", 
                                         names(annotated.utr3)[1]))
            utr3end <- 
                as.numeric(gsub("^.*_SEP_", "", 
                                names(annotated.utr3)[length(annotated.utr3)]))
            last.annotated.utr3 <- annotated.utr3[length(annotated.utr3)]
            ## drop low coverage from the end
            if(length(next.exon.gap)>1){
                next.exon.gap.abun <- 
                    cumsum(rev(next.exon.gap))/1:length(next.exon.gap)
                id <- which(next.exon.gap.abun>long_coverage_threshold)
                if(length(id)>0){
                    next.exon.gap.abun.1 <- length(next.exon.gap) - id[1]
                }else{
                    next.exon.gap.abun.1 <- 0
                }   
            }else{
                next.exon.gap.abun.1 <- 0
            }
            ## search for polyA site
            if(search_distal_polyA_end &&
                   next.exon.gap.abun.1>100 && 
                   class(classifier)=="PASclassifier") {
                next.exon.gap.abun.1 <- 
                    searchPolyAfromEnd(classifier, classifier_cutoff, 
                                       next.exon.gap[1:next.exon.gap.abun.1], 
                                       as.character(seqnames(curr_UTR.ele))[1],
                                       as.character(strand(curr_UTR.ele))[1], 
                                       genome, window_size)
            }
            
            saved.proximal.apa <- 
                as.numeric(gsub("^.*_SEP_", "", 
                                names(annotated.utr3)[length(annotated.utr3)]))
            saved.id <- length(annotated.utr3)
            if(next.exon.gap.abun.1 > window_size){
                type <- "novalDistal"
                Predicted_Distal_APA <- 
                    as.numeric(gsub("^.*_SEP_", "", 
                                    names(next.exon.gap)[
                                        next.exon.gap.abun.1]))
                .ele <- 
                    rbind(chr.cov.merge.ele[grepl("utr3", 
                                                  rownames(chr.cov.merge.ele)), 
                                            , drop=FALSE],
                              chr.cov.merge.ele[
                                  names(next.exon.gap)[1:next.exon.gap.abun.1], 
                                  , drop=FALSE])
                
            }else{
                type <- "distal"
                Predicted_Distal_APA <- 
                    as.numeric(gsub("^.*_SEP_", "", 
                                    names(annotated.utr3)[
                                        length(annotated.utr3)]))
                .ele <- 
                    chr.cov.merge.ele[grepl("utr3", 
                                            rownames(chr.cov.merge.ele)), 
                                      , drop=FALSE]
            }
            if(grepl("proximalCP", 
                     curr_UTR.ele[curr_UTR.ele$id=="utr3"]$feature[1])){
                ##check proximal CP site position is in the range or not
                Predicted_Proximal_APA <- 
                    as.integer(
                        unlist(
                        strsplit(
                        curr_UTR.ele[curr_UTR.ele$id=="utr3"]$feature[1], "_")
                        )[-1])[1]
                fit_value <- NA
            }else{
                if(nrow(.ele) > max(c(search_point_START, MINSIZE))){
                    if(!is.na(cutEnd)){
                        if(cutEnd<1){
                            .ele <- 
                                .ele[1:floor((nrow(.ele)-1)*(1-cutEnd)), 
                                     , drop=FALSE]
                        }else{
                            .ele <- 
                                .ele[1:(length(.ele)-1-floor(cutEnd)), 
                                     , drop=FALSE]
                        }
                    }
                    
                    rownames(.ele) <- gsub("^.*_SEP_", "", rownames(.ele))
                    .l <- nrow(.ele)
                    if(is.na(search_point_END)) {
                        search_point_end <- .l - 1
                    }else{
                        if(search_point_END<0) {
                            search_point_end <- floor(.l + search_point_END)
                        }else{
                            if(search_point_END<1) {
                                search_point_end<- 
                                    ceiling(.l*(1-search_point_END))
                            }else {
                                search_point_end <- 
                                    floor(.l - search_point_END)
                            }
                        }
                    }
                    if(search_point_end > search_point_START){
                        fos <- apply(.ele, 2, optimalSegmentation, 
                                     search_point_START=search_point_START, 
                                     search_point_END=search_point_end)
                        cov_diff <- sapply(fos, "[[", "cov_diff")
                        cov_diff <- rowMeans(cov_diff)
                        
                        ## select the best 1
                        idx <- valley(cov_diff, search_point_START, 
                                      search_point_end, n=-1, savedID=saved.id)
                        if(search_point_START<MINSIZE) 
                            idx <- idx[idx!=search_point_START]
                        idx1 <- if(length(idx)>0) idx[1] else MINSIZE
                        
                        seqname <- as.character(seqnames(curr_UTR.ele))[1]
                        strand <- as.character(strand(curr_UTR.ele))[1]
                        start <- ifelse(strand=="+", 
                                        as.numeric(rownames(.ele)[1]), 
                                        as.numeric(rownames(.ele)[nrow(.ele)]))
                        idx <- filterIdx(idx, idx1, cov_diff, PolyA_PWM, 
                                         seqname, start, strand, genome, 
                                         classifier, classifier_cutoff, 
                                         shift_range, search_point_START, 
                                         type, saved.id, window_size=NA)
                        
                        fit_value <- cov_diff[idx]
                        Predicted_Proximal_APA <- 
                            as.numeric(rownames(.ele)[idx])
                    }else{
                        fit_value <- NA
                        Predicted_Proximal_APA <- NA
                    }
                }else{
                    fit_value <- NA
                    Predicted_Proximal_APA <- NA
                } 
            }                   
            
            list(fit_value=fit_value,
                 Predicted_Proximal_APA=Predicted_Proximal_APA,
                 Predicted_Distal_APA=Predicted_Distal_APA,
                 utr3start=utr3start,
                 utr3end=utr3end,
                 type=type)
        }, chr.cov.merge, 
        conn_next_utr3, 
        curr_UTR[names(chr.cov.merge)], 
        SIMPLIFY=FALSE)
        chr.abun <- do.call(rbind, chr.abun)
    }else{
        NULL
    }
}

filterIdx <- function(idx, idx1, cov_diff, PolyA_PWM, 
                      seqname, start, strand, genome, 
                      classifier, classifier_cutoff, 
                      shift_range, search_point_START, 
                      type, saved.id, window_size=NA){
    if(class(PolyA_PWM)=="matrix" && length(idx)>0){
        if(is.na(window_size)){
            pos <- if(strand=="+") start + idx - 1 else start - idx + 1
        }else{
            pos <- 
                if(strand=="+"){
                    start + floor((idx - .5) * window_size)
                }else{
                    start - floor((idx - .5) * window_size)
                }
        }
        idx <- PAscore(seqname, pos, str=strand, 
                       idx, PWM=PolyA_PWM, genome=genome) 
    }
    if(class(classifier)=="PASclassifier" && length(idx)>0){
        if(is.na(window_size)){
            if(shift_range>10){
                idx_lo <- idx - shift_range
                idx_up <- idx + shift_range
                idx <- as.integer(mapply(function(a, b, c) 
                    unique(sort(c(seq(a, b, by=10), c))), idx_lo, idx_up, idx))
                idx <- 
                    unique(idx[idx>=search_point_START & 
                                   idx < length(cov_diff)])
            }
            pos <- if(strand=="+") start + idx - 1 else start - idx + 1
        }else{
            pos <- if(strand=="+"){ 
                    start + floor((idx - .5) * window_size) 
                }else{
                    start - floor((idx - .5) * window_size)
                }
        }
        idx <- PAscore2(seqname, pos, strand, idx, 
                        genome, classifier, classifier_cutoff)
    }
    
    if(length(idx)==0){
        if(type=="novalDistal"){
            idx <- saved.id
        }else{
            idx <- idx1
        }
    }else{
        ## re-sort the idx
        ## cov_diff[idx] lower is better (already filtered by valley)
        ## ratio of two fragment smaller is letter
        ## ratio <- abs(log2(idx) - log2(length(cov_diff)-idx))
        ## ratio <- cut(ratio, breaks=c(0, 3.321928, 6.643856, 9.965784, Inf), 
        ##             labels=c(1, 2, 3, 4))
        ##ratio <- as.numeric(as.character(ratio))
        ##idx <- idx[order(ratio)]
        idx <- idx[1]
    }
    idx
}

optimalSegmentation <- function(.ele, search_point_START, 
                                search_point_END, n=1, savedID=NA){
    .l <- length(.ele)
    short_UTR_abun <- cumsum(.ele)/1:.l
    long_UTR_abun <- cumsum(rev(.ele))/1:.l
    long_UTR_abun <- rev(long_UTR_abun)
    short_UTR_abun <- short_UTR_abun[-length(short_UTR_abun)]
    long_UTR_abun <- long_UTR_abun[-1]
    short_UTR_abun <- short_UTR_abun - long_UTR_abun
    short_UTR_abun <- ifelse(short_UTR_abun<0, 0, short_UTR_abun)
    cov_diff <- numeric(length(short_UTR_abun))
    ss <- max(search_point_START, 1)
    se <- min(search_point_END, .l)
    for(i in ss:se){
        cov_diff_tmp <- .ele
        cov_diff_tmp <- cov_diff_tmp-long_UTR_abun[i]
        cov_diff_tmp[1:i] <- cov_diff_tmp[1:i] - short_UTR_abun[i]
        cov_diff[i] <- mean(cov_diff_tmp^2)
    }
    idx <- valley(cov_diff, ss, se, n, savedID)
    return(list(cov_diff=cov_diff, idx=idx))
}

PAscore <- function(seqname, pos, str, idx, PWM, genome, ups=50, dws=50){
    if(length(pos)<1){
        return(NULL)
    }
    start <- pos-ups
    start[start<1] <- 1
    end <- pos + dws
    gr <- GRanges(seqname, IRanges(start, end, names=as.character(pos)), 
                  strand=str)
    seq <- getSeq(genome, gr)
    mT <- lapply(seq, matchPWM, pwm=PWM, min.score="70%", with.score=TRUE)
    hits <- sapply(mT, function(.ele) {
        if(class(.ele)!="XStringViews") return(FALSE)
        if(length(.ele)==0) return(FALSE)
        TRUE
    })
    idx[hits]
}


PAscore2 <- function(seqname, pos, str, idx, 
                     genome, classifier, classifier_cutoff){
    if(length(pos)<1){
        return(NULL)
    }
    gr <- 
        GRanges(seqname, 
                IRanges(pos, pos, names=as.character(pos)), 
                strand=str)
    testSet.NaiveBayes <- 
        buildFeatureVector(gr, BSgenomeName = genome, 
                           upstream = classifier@info@upstream,
                           downstream = classifier@info@downstream, 
                           wordSize = classifier@info@wordSize, 
                           alphabet=classifier@info@alphabet,
                           sampleType = "unknown",replaceNAdistance = 30, 
                           method = "NaiveBayes", ZeroBasedIndex = 1, 
                           fetchSeq = TRUE)
    suppressMessages(pred.prob.test <- 
                         predictTestSet(testSet.NaiveBayes=testSet.NaiveBayes, 
                                        classifier=classifier, 
                                        outputFile=NULL, 
                                        assignmentCutoff=classifier_cutoff))
    idx[pred.prob.test[, "pred.class"]==1]
}

##remove utr3---___---utr3, need to improve.
removeUTR3__UTR3 <- function(x){
    ## smooth by window_size=10
    ws <- 10
    len <- length(x)
    if(len>100){
        y <- rowsum(x, group=rep(1:ceiling(len/ws), each=ws)[1:len], 
                    reorder=FALSE)
        y <- y[-length(y)]
        id <- valley(y, 1, length(y), 1, filterByPval=FALSE)
        if(length(id)>0){
            if(id<3) id <- 3
            pval <- try(t.test(y[1:(id-2)], y[-(1:(id+2))])$p.value, 
                        silent=TRUE)
            if(class(pval)=="numeric" && length(pval)==1 && !is.na(pval)){
                if(pval < 0.001){
                    x <- x[1:((id-1)*ws)]
                }
            }
        }
    }
    x
}

searchPolyAfromEnd <- function(classifier, classifier_cutoff, 
                               gap.cov, chrom, strand, genome, 
                               window_size){
    coor <- as.integer(gsub("^.*_SEP_", "", names(gap.cov)))
    start <- coor[length(coor)]
    end <- ifelse(length(coor) > 2*window_size, 
                  coor[length(coor)-2*window_size], 
                  coor[1])
    pos <- seq(start, end, by=ifelse(strand=="+", -10, 10))
    idx <- match(pos, coor)
    idx1 <- PAscore2(chrom, pos, strand, idx, genome, 
                     classifier, classifier_cutoff)
    if(length(idx1)>0) {
        idx1 <- idx1[1]
    }else{
        idx1 <- idx[1]
    }
    idx1
}