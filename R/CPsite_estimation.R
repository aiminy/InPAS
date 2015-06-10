CPsite_estimation <- function(chr.cov, utr3, MINSIZE, window_size, 
                              search_point_START, search_point_END, 
                              cutStart, cutEnd, adjust_distal_polyA_end,  
                              background, z2s,
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
                                    long_coverage_threshold,
                                    background, 
                                    z2s)
        ##step2 adjust distal cp sites
        if(adjust_distal_polyA_end && class(classifier)=="PASclassifier") 
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