searchDistalCPs <- function(chr.cov.merge,
                            conn_next_utr3,
                            curr_UTR, 
                            window_size,
                            depth.weight,
                            long_coverage_threshold,
                            background,
                            z2s){
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
        ## stop and any connected two windows different from 30 times, 
        ## to separate connected two UTRs or strange peaks.
        if(length(next.exon.gap)>1){
            next.exon.gap.ids <- 
                rep(1:ceiling(length(next.exon.gap)/window_size), 
                    each=window_size)[1:length(next.exon.gap)]
            next.exon.gap.split <- split(next.exon.gap, next.exon.gap.ids)
            next.exon.gap.split <- sapply(next.exon.gap.split, mean)
            next.exon.gap.split.diff <- 
                log2(next.exon.gap.split[-1]+0.01) - 
                log2(next.exon.gap.split[-length(next.exon.gap.split)]+0.01)
            id <- which(abs(next.exon.gap.split.diff)>log2(30))
            ### drop by background by window
            ## the gap == 2Xwindow size low coverage
            cont.num <- function(x){
                if(length(x)==1) return(x)
                r <- which(diff(x)==1)
                if(length(r)>0){
                    return(x[r[1]+1])
                }else{
                    return(x[length(x)])
                }
            }
            if(background!="same_as_long_coverage_threshold"){
                z2_threshold <- z2s[names(curr_UTR.ele[curr_UTR.ele$id=="utr3"])]
                if(!is.na(z2_threshold)){
                    z2 <- which(next.exon.gap.split < z2_threshold)
                    if(length(z2)>0){
                        z2 <- cont.num(z2)
                        id <- min(c(z2, id))
                    }
                }
            }else{
                ##background == long_coverage_threshold
                z2 <- which(next.exon.gap.split < long_coverage_threshold)
                if(length(z2)>0){
                    z2 <- cont.num(z2)
                    id <- min(c(z2, id))
                }
            }
            if(length(id)>0){
                conn_next_utr <- FALSE
                id <- id[1] - 1
                if(id>0){
                    next.exon.gap <- next.exon.gap[1:(window_size*id)]
                }else{
                    next.exon.gap <- numeric(0)
                }
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
        ### drop by long_coverage_threshold
        if(length(next.exon.gap)>5){
            ## make it simple, smooth and cutoff
            fit.loess <- loess.smooth(1:length(next.exon.gap), 
                                      next.exon.gap, 
                                      evaluation=length(next.exon.gap),
                                      family='gaussian', 
                                      span=.75, degree=1)
            next.exon.gap.abun <- fit.loess$y
            id <- which(next.exon.gap.abun>=long_coverage_threshold)
            if(length(id)>0){
                next.exon.gap.abun <- next.exon.gap[1:id[length(id)]]
                id <- which(next.exon.gap.abun>=long_coverage_threshold)
            }
            if(length(id)>0){
                next.exon.gap.abun.1 <- id[length(id)]
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