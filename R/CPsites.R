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
        if(length(sn)<=p) return(sn)
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
                    cutStart=window_size, cutEnd=0, search_distal_polyA_end=FALSE, 
                    coverage_threshold=5, long_coverage_threshold=2, 
                    gcCompensation=NA, mappabilityCompensation=NA, 
                    FFT=FALSE, fft.sm.power=20, 
                    PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, 
                    shift_range=window_size, BPPARAM=NULL){
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
                unname(shorten_UTR_estimation[!sapply(shorten_UTR_estimation, 
                                                      is.null)]))
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
            if(length(dim(.covList))==0) .covList <- t(.covList)
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
        ##step1 search distal cp sites
        curr_UTR <- curr_UTR[names(chr.cov.merge)]
        chr.abun <- searchDistalCPs(chr.cov.merge,
                                    conn_next_utr3,
                                    curr_UTR,
                                    window_size,
                                    depth.weight,
                                    long_coverage_threshold)
        ##step2 adjust distal cp sites
        if(search_distal_polyA_end && class(classifier)=="PASclassifier") 
            chr.abun <- distalAdj(chr.abun, classifier, classifier_cutoff,
                                  shift_range, genome)
        
        ##step3 search proximal cp sites
        chr.abun <- searchProximalCPs(chr.abun, curr_UTR, 
                                      window_size, MINSIZE,
                                      cutEnd, 
                                      search_point_START, search_point_END)
        
        
        ##step4 adjust proximal cp sites
        if(class(PolyA_PWM)=="matrix" || class(classifier)=="PASclassifier"){
            chr.abun <- proximalAdj(chr.abun, MINSIZE, PolyA_PWM, 
                                    genome, classifier, classifier_cutoff, 
                                    shift_range, search_point_START)
        }
        
        chr.abun <- polishCPs(chr.abun)
    }else{
        NULL
    }
}

## sub functions
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

searchDistalCPs <- function(chr.cov.merge,
                            conn_next_utr3,
                            curr_UTR, 
                            window_size,
                            depth.weight,
                            long_coverage_threshold){
    distalCPs <- mapply(function(chr.cov.merge.ele, 
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
        
        info <- as.data.frame(curr_UTR.ele)[1, ]
        info$utr3start <- utr3start
        info$utr3end <- utr3end
        info$distalCP <- next.exon.gap.abun.1
        list(info=info, 
             cov=chr.cov.merge.ele,
             gap=next.exon.gap,
             annotated.utr3=annotated.utr3
        )
    }, chr.cov.merge, 
    conn_next_utr3, 
    curr_UTR, 
    SIMPLIFY=FALSE)
    dCPs <- do.call(rbind, lapply(distalCPs, `[[`, "info"))
    chr.cov.merge <- lapply(distalCPs, `[[`, "cov")
    next.exon.gap <- lapply(distalCPs, `[[`, "gap")
    annotated.utr3 <- lapply(distalCPs, `[[`, "annotated.utr3")
    list(dCPs=dCPs, chr.cov.merge=chr.cov.merge,
         next.exon.gap=next.exon.gap,
         annotated.utr3=annotated.utr3)
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

PAscore2 <- function(seqname, pos, str, idx, idx.gp,
                     genome, classifier, classifier_cutoff){
    if(length(pos)<1){
        return(NULL)
    }
    coor <- paste(seqname, pos, str, sep="_")
    gr <- 
        GRanges(seqname, 
                IRanges(pos, pos, names=coor), 
                strand=str)
    gr$id <- 1:length(gr)
    coor.id <- !duplicated(coor)
    gr$duplicated <- gr$id[match(coor, coor[coor.id])]
    gr.s <- gr[coor.id]
    testSet.NaiveBayes <- 
        buildFeatureVector(gr.s, BSgenomeName = genome, 
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
    if(any(duplicated(coor))){
        ## need to recover the order of inputs
        pred.prob.test <- pred.prob.test[gr$duplicated, , drop=FALSE]
        pred.prob.test[, "PeakName"] <- names(gr)
    }
    pred.prob.test <- cbind(pred.prob.test[,1:4], idx, idx.gp)
    pred.prob.test <- pred.prob.test[pred.prob.test[, "pred.class"]==1, ]
    pred.prob.test <- pred.prob.test[order(pred.prob.test[, "idx.gp"],
                                           -pred.prob.test[, "prob True"]), ]
    pred.prob.test <- pred.prob.test[!duplicated(pred.prob.test[, "idx.gp"]), ]
    pred.prob.test
}

distalAdj <- function(distalCPs, classifier, classifier_cutoff,
                      shift_range, genome){
    dCPs <- distalCPs$dCPs
    next.exon.gap <- distalCPs$next.exon.gap
    gap.cov <- mapply(function(gap, cp, ID, strand){
        if(cp>0){
            coor <- as.integer(gsub("^.*_SEP_", "", names(gap[1:cp])))
            start <- coor[length(coor)]
            end <- ifelse(length(coor) > 2*shift_range, 
                          coor[length(coor)-2*shift_range], 
                          coor[1])
            pos <- seq(start, end, by=ifelse(strand=="+", -10, 10))
            idx <- match(pos, coor)
            cbind(pos, idx, ID)
        }else{
            NULL
        }
    }, next.exon.gap, dCPs$distalCP, 
    1:length(next.exon.gap), dCPs$strand,
    SIMPLIFY=FALSE)
    gap.cov <- do.call(rbind, gap.cov)
    if(length(gap.cov)>0){
        idx <- PAscore2(dCPs$seqnames[gap.cov[, "ID"]],
                        gap.cov[, "pos"],
                        dCPs$strand[gap.cov[, "ID"]],
                        gap.cov[, "idx"],
                        gap.cov[, "ID"],
                        genome, classifier, classifier_cutoff)
        distalCPs$dCPs[idx$idx.gp, "distalCP"] <- idx$idx
    }
    distalCPs
}

searchProximalCPs <- function(CPs, curr_UTR, 
                              window_size, MINSIZE,
                              cutEnd, 
                              search_point_START, search_point_END){
    dCPs <- CPs$dCPs
    chr.cov.merge <- CPs$chr.cov.merge
    next.exon.gap <- CPs$next.exon.gap
    annotated.utr3 <- CPs$annotated.utr3
    saved.id <- dCPs$length <- sapply(annotated.utr3, length)
    saved.proximal.apa <- mapply(function(.ele, .len){
        as.numeric(gsub("^.*_SEP_", "", names(.ele)[.len]))
    }, annotated.utr3, saved.id)
    flag <- dCPs$distalCP > window_size
    dCPs$length[flag] <- dCPs$length[flag] + dCPs$distalCP[flag]
    dCPs$type <- ifelse(flag, "novalDistal", "distal")
    dist_apa <- function(d, id){
        ifelse(id>0, 
               as.numeric(rownames(d)[id]), 
               0)
    }
    chr.cov.merge <- lapply(chr.cov.merge, function(.ele){
        rownames(.ele) <- gsub("^.*_SEP_", "", rownames(.ele))
        .ele
    })
    dCPs$Predicted_Distal_APA <- mapply(dist_apa, chr.cov.merge, dCPs$length)
    chr.cov.merge <- mapply(function(d, len){
        if(len>0) {
            d[1:len, , drop=FALSE]
        }else{
            d[FALSE,, drop=FALSE]
        }
    }, chr.cov.merge, dCPs$length, SIMPLIFY=FALSE)
    proximalCP <- sapply(curr_UTR, function(.ele) 
        grepl("proximalCP", .ele[.ele$id=="utr3"]$feature[1]))
    Predicted_Proximal_APA <- vector("list", length=nrow(dCPs))
    fit_value <- vector("list", length=nrow(dCPs))
    dCPs$fit_value <- NA
    if(sum(proximalCP)>0){
        Predicted_Proximal_APA[proximalCP] <- 
            lapply(curr_UTR[proximalCP], function(.ele){
                as.integer(unlist(strsplit(
                    .ele[.ele$id=="utr3"]$feature[1], "_")[2]))
            })
    }
    if(!is.na(cutEnd)){
        if(cutEnd<1){
            chr.cov.merge <- lapply(chr.cov.merge, function(.ele){
                .ele[1:floor((nrow(.ele)-1)*(1-cutEnd)), 
                     , drop=FALSE]
            })
        }else{
            chr.cov.merge <- lapply(chr.cov.merge, function(.ele){
                .ele[1:(length(.ele)-1-floor(cutEnd)), 
                     , drop=FALSE]
            })
        }
    }
    
    minStartPos <- dCPs$length >= max(c(search_point_START, MINSIZE))
    len <- sapply(chr.cov.merge, nrow)
    search_point_END <- rep(abs(search_point_END), nrow(dCPs))
    search_point_end <- ifelse(is.na(search_point_END),
                               len - 1,
                               ifelse(search_point_END<1, 
                                      ceiling(len*(1-search_point_END)),
                                      floor(len - search_point_END)))
    flag <- minStartPos & 
        (search_point_end > search_point_START) & 
        (!proximalCP)
    fit_value[flag] <- mapply(function(.ele, search_point_END){
        fos <- apply(.ele, 2, optimalSegmentation, 
                     search_point_START=search_point_START, 
                     search_point_END=search_point_END)
        cov_diff <- sapply(fos, "[[", "cov_diff")
        cov_diff <- rowMeans(cov_diff)
    }, chr.cov.merge[flag], search_point_end[flag], SIMPLIFY=FALSE)
    Predicted_Proximal_APA[flag] <- 
        mapply(function(cov_diff, search_point_END, savedID){
            idx <- valley(cov_diff, search_point_START, 
                          search_point_END, n=-1, savedID=savedID)
            if(search_point_START<MINSIZE) 
                idx <- idx[idx!=search_point_START]
            idx
        }, fit_value[flag], search_point_end[flag], saved.id[flag], 
        SIMPLIFY=FALSE)
    
    
    idx1 <- lapply(Predicted_Proximal_APA, `[`, 1)
    idx1[sapply(idx1, length)==0] <- NA
    idx1 <- unlist(idx1)
    
    list(dCPs=dCPs, chr.cov.merge=chr.cov.merge,
         next.exon.gap=next.exon.gap,
         annotated.utr3=annotated.utr3,
         flag=flag, fit_value=fit_value, 
         Predicted_Proximal_APA=Predicted_Proximal_APA,
         saved.id=saved.id, idx1=idx1)
}

proximalAdjByPWM <- function(idx, PolyA_PWM, seqnames, starts, strands,
                             genome, shift_range, search_point_START){
    mapply(function(id, seqname, start, strand){
        if(length(id)>0){
            pos <- if(strand=="+") start + id - 1 else start - id + 1
            id <- PAscore(seqname, pos, strand, 
                          id, PWM=PolyA_PWM, genome=genome,
                          ups=shift_range+25,
                          dws=shift_range+25)
        }
        id
    }, idx, seqnames, starts, strands, SIMPLIFY=FALSE)
}

proximalAdjByCleanUpdTSeq <- function(idx.list, cov_diff.list, 
                                      seqnames, starts, strands, 
                                      genome, classifier, classifier_cutoff,
                                      shift_range, search_point_START){
    idx.len <- sapply(idx.list, length)
    offsite <- 10^nchar(as.character(max(idx.len)))
    pos.matrix <- mapply(function(idx, start, strand, cov_diff, ID){
        if(length(idx)==0) return(NULL)
        if(is.na(idx[1])) return(NULL)
        idx.gp <- 1:length(idx)
        if(shift_range>10){
            idx_lo <- idx - shift_range
            idx_up <- idx + shift_range
            idx <- mapply(function(a, b, c) 
                unique(sort(c(seq(a, b, by=10), c))), 
                idx_lo, idx_up, idx, SIMPLIFY=FALSE)
            idx.gp <- rep(1:length(idx_lo), sapply(idx, length))
            idx <- as.integer(unlist(idx))
            idx <- cbind(idx, idx.gp)
            idx <- idx[idx[,1]>=search_point_START & 
                           idx[,1] < length(cov_diff), , drop=FALSE]
            idx <- idx[!duplicated(idx[,1]), , drop=FALSE]
            idx.gp <- idx[, 2]
            idx <- idx[, 1]
        }
        if(length(idx)>0){
            pos <- if(strand=="+") start + idx - 1 else start - idx + 1
            cbind(pos, idx, idx.gp, ID)
        }else{
            NULL
        }
    }, idx.list, starts, strands, cov_diff.list, 
    1:length(idx.list), SIMPLIFY=FALSE)
    pos.matrix <- do.call(rbind, pos.matrix)
    
    if(length(pos.matrix)>0){
        idx <- PAscore2(seqnames[pos.matrix[, "ID"]],
                        pos.matrix[, "pos"],
                        strands[pos.matrix[, "ID"]],
                        pos.matrix[, "idx"],
                        pos.matrix[, "ID"]*offsite + pos.matrix[, "idx.gp"],
                        genome, classifier, classifier_cutoff)
        idx$ID <- floor(idx$idx.gp/offsite)
        idx <- idx[!duplicated(idx$ID), ]
        IDs <- unique(pos.matrix[, "ID"])
        idx.idx <- idx$idx[match(IDs, idx$ID)]
        idx.idx[is.na(idx.idx)] <- sapply(idx.list[IDs[is.na(idx.idx)]], `[`, 1)
        idx.list[IDs] <- idx.idx
    }
    idx.list
}

proximalAdj <- function(CPs, MINSIZE, PolyA_PWM, 
                        genome, classifier, classifier_cutoff, 
                        shift_range, search_point_START){
    dCPs <- CPs$dCPs
    flag <- CPs$flag
    seqnames <- as.character(dCPs$seqnames)
    strands <- as.character(dCPs$strand)
    starts <- coors <- 
        lapply(CPs$chr.cov.merge, function(.ele) as.numeric(rownames(.ele)))
    starts[strands=="-"] <- lapply(starts[strands=="-"], rev)
    starts <- sapply(starts, `[`, 1)
    idx.list <- CPs$Predicted_Proximal_APA
    if(class(PolyA_PWM)=="matrix"){
        idx.list <- proximalAdjByPWM(idx.list, PolyA_PWM, seqnames, starts,
                                     strands, genome, shift_range, 
                                     search_point_START)
    }
    cov_diff.list <- CPs$fit_value
    if(class(classifier)=="PASclassifier"){
        idx.list <- proximalAdjByCleanUpdTSeq(idx.list, cov_diff.list, 
                                              seqnames, starts, strands, 
                                              genome, 
                                              classifier, classifier_cutoff,
                                              shift_range, search_point_START)
    }
    CPs$Predicted_Proximal_APA[flag] <- idx.list[flag]
    
    CPs
}

polishCPs <- function(CPs){
    proximal.apa.len <- sapply(CPs$Predicted_Proximal_APA, length)
    CPs$Predicted_Proximal_APA[proximal.apa.len==0] <- 
        ifelse(CPs$dCPs$type[proximal.apa.len==0]=="distal",
               CPs$idx1[proximal.apa.len==0],
               CPs$saved.id[proximal.apa.len==0])
    coors <- lapply(CPs$chr.cov.merge, function(.ele) as.numeric(rownames(.ele)))
    flag <- CPs$flag
    CPs$fit_value[flag] <- mapply(function(cov_diff, idx){
        ifelse(is.na(idx[1]) | is.null(idx[1]), NA, cov_diff[idx[1]])
    }, CPs$fit_value[flag], 
    CPs$Predicted_Proximal_APA[flag], 
    SIMPLIFY=FALSE)
    CPs$Predicted_Proximal_APA <- mapply(function(.ele, .id){
        ifelse(is.na(.id[1]) | is.null(.id[1]), NA, .ele[.id[1]])
    }, coors,
    CPs$Predicted_Proximal_APA,
    SIMPLIFY=FALSE)
    CPs$fit_value[sapply(CPs$fit_value, length)==0] <- NA
    CPs$dCPs$fit_value <- unlist(CPs$fit_value)
    CPs$Predicted_Proximal_APA[sapply(CPs$Predicted_Proximal_APA, length)==0] <- NA
    CPs$dCPs$Predicted_Proximal_APA <- unlist(CPs$Predicted_Proximal_APA)
    CPs$dCPs[, c("fit_value", "Predicted_Proximal_APA", 
                 "Predicted_Distal_APA", "utr3start", 
                 "utr3end", "type")]
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