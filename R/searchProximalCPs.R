searchProximalCPs <- function(CPs, curr_UTR, 
                              window_size, MINSIZE,
                              cutEnd, 
                              search_point_START, 
                              search_point_END,
                              two_way=FALSE){
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
    dCPs$type <- ifelse(flag, "novel distal", "novel proximal")
    dCPs$type[grepl("proximalCP", dCPs$annotatedProximalCP) & !flag] <- 
        "annotated proximal"
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
        grepl("proximalCP", .ele[.ele$feature=="utr3"]$annotatedProximalCP[1]))
    Predicted_Proximal_APA <- vector("list", length=nrow(dCPs))
    fit_value <- vector("list", length=nrow(dCPs))
    Predicted_Proximal_APA_rev <- vector("list", length=nrow(dCPs))
    fit_value_rev <- vector("list", length=nrow(dCPs))
    dCPs$fit_value <- NA
    if(sum(proximalCP)>0){
        Predicted_Proximal_APA[proximalCP] <- 
            lapply(curr_UTR[proximalCP], function(.ele){
                as.integer(unlist(strsplit(
                    .ele[.ele$feature=="utr3"]$annotatedProximalCP[1], "_")[2]))
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
                .ele[1:max(nrow(.ele)-1-floor(cutEnd), 1), 
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
    ##forward
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
                          search_point_END, n=5, savedID=savedID)
            if(search_point_START<MINSIZE) 
                idx <- idx[idx!=search_point_START]
            idx
        }, fit_value[flag], search_point_end[flag], saved.id[flag], 
        SIMPLIFY=FALSE)
    if(two_way){
        ##reverse
        fit_value_rev[flag] <- mapply(function(.ele, search_point_END){
            nr <- nrow(.ele)
            fos <- apply(.ele[nr:1,,drop=FALSE], 2, optimalSegmentation, 
                         search_point_START=nr-search_point_END, 
                         search_point_END=nr-search_point_START)
            cov_diff <- sapply(fos, "[[", "cov_diff")
            cov_diff <- rowMeans(cov_diff)
            cov_diff <- cov_diff[length(cov_diff):1]
        }, chr.cov.merge[flag], search_point_end[flag], SIMPLIFY=FALSE)
        Predicted_Proximal_APA_rev[flag] <- 
            mapply(function(cov_diff, search_point_END, savedID){
                idx <- valley(cov_diff, search_point_START, 
                              search_point_END, n=5, savedID=savedID)
                if(search_point_START<MINSIZE) 
                    idx <- idx[idx!=search_point_START]
                idx
            }, fit_value_rev[flag], search_point_end[flag], saved.id[flag], 
            SIMPLIFY=FALSE)
        ##combine forward and reverse
        Predicted_Proximal_APA[flag] <- mapply(function(fv, idx, fv_rev, idx_rev){
            ## sort the results
            f <- fv[idx]
            names(f) <- idx
            idx_rev <- idx_rev[!idx_rev %in% idx]
            r <- fv_rev[idx_rev]
            names(r) <- idx_rev
            fr <- sort(c(f, r), decreasing = FALSE)
            as.numeric(names(fr))[1:min(length(fr), 6)]
        }, fit_value[flag], Predicted_Proximal_APA[flag], 
        fit_value_rev[flag], Predicted_Proximal_APA_rev[flag], 
        SIMPLIFY = FALSE)
    }
    
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
