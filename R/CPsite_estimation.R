CPsite_estimation <- function(chr.cov, utr3, MINSIZE, window_size, 
                              search_point_START, search_point_END, 
                              cutStart, cutEnd, adjust_distal_polyA_end,  
                              background, z2s,
                              coverage_threshold, long_coverage_threshold, 
                              PolyA_PWM, classifier, classifier_cutoff, 
                              shift_range, depth.weight, genome, 
                              tmpfolder=NULL, silence=TRUE){
    chr.cov <- chr.cov[sapply(chr.cov, mean)>0]
    if(length(chr.cov)==0){
        return(NULL)
    }
    curr_UTR.gr <- utr3[names(chr.cov)]
    seqn <- as.character(seqnames(curr_UTR.gr))[1]
    
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
        chr.utr3TotalCov <- 
            mapply(function(.covList, .start, .end, .property){
                #set names for each position 
                .posList <- .start:.end
                if(length(dim(.covList))==0) .covList <- t(.covList)
                rownames(.covList) <- 
                    paste(.property, .posList, sep="_SEP_")
                .covList
            }, chr.utr3TotalCov, start(.UTR), end(.UTR), 
            .UTR$id, SIMPLIFY=FALSE)
        chr.utr3TotalCov <- do.call(rbind, chr.utr3TotalCov)
        if(as.character(strand(.UTR))[1] == "-"){##reverse the negative strand
            chr.utr3TotalCov <- 
                chr.utr3TotalCov[rev(rownames(chr.utr3TotalCov)), 
                                 , drop=FALSE]
        }
        if(!is.na(cutStart)){
            if(cutStart<1) 
                cutStart <- floor(length(chr.utr3TotalCov)*cutStart)
            if(cutStart>0) 
                chr.utr3TotalCov <- 
                chr.utr3TotalCov[-(1:cutStart), , drop=FALSE]
        }
        
        chr.utr3TotalCov
    })
    if(!silence) message("chromsome", seqn, "coverage merged.")
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
    if(!silence) message("chromsome", seqn, "quality filtered by first 100nt.")
    
    if(!is.null(tmpfolder)){
        if(file.exists(file.path(tmpfolder, seqn))){
            load(file.path(tmpfolder, seqn))
            if(!silence) message("chromsome", seqn, "coverage loaded.")
        }
    }
    if(!exists("CPsite_search_step")){
        CPsite_search_step <- 0
    }
    
    if(length(chr.cov.merge)>0){
        ##step1 search distal cp sites
        if(!silence) message("chromsome", seqn, "distal search ... start")
        curr_UTR <- curr_UTR[names(chr.cov.merge)]
        if(CPsite_search_step<1){
            CPsite_search_step <- 1
            chr.abun <- searchDistalCPs(chr.cov.merge,
                                        conn_next_utr3,
                                        curr_UTR,
                                        window_size,
                                        depth.weight,
                                        long_coverage_threshold,
                                        background, 
                                        z2s)
            if(!is.null(tmpfolder)) 
                save(list=c("chr.abun", "CPsite_search_step"),
                     file=file.path(tmpfolder, seqn))
            if(!silence) message("chromsome", seqn, "distal search ... done")
        }
        
        ##step2 adjust distal cp sites
        if(!silence) message("chromsome", seqn, "distal adjust ... start")
        if(adjust_distal_polyA_end && class(classifier)=="PASclassifier"){
            if(CPsite_search_step < 2){
                CPsite_search_step <- 2
                chr.abun <- distalAdj(chr.abun, classifier, classifier_cutoff,
                                      shift_range, genome)
                if(!is.null(tmpfolder)) 
                    save(list=c("chr.abun", "CPsite_search_step"),
                         file=file.path(tmpfolder, seqn))
                if(!silence) 
                    message("chromsome", seqn, "distal adjust ... done")
            }
        }  
        
        ##step3 search proximal cp sites
        if(!silence) message("chromsome", seqn, "proximal search ... start")
        if(CPsite_search_step < 3){
            CPsite_search_step <- 3
            chr.abun <- searchProximalCPs(chr.abun, curr_UTR, 
                                          window_size, MINSIZE,
                                          cutEnd, 
                                          search_point_START, search_point_END)
            if(!is.null(tmpfolder)) 
                save(list=c("chr.abun", "CPsite_search_step"),
                     file=file.path(tmpfolder, seqn))
            if(!silence) 
                message("chromsome", seqn, "proximal searched ... done")
        }
        
        ##step4 adjust proximal cp sites
        if(!silence) message("chromsome", seqn, "proximal adjust ... start")
        if(class(PolyA_PWM)=="matrix" || class(classifier)=="PASclassifier"){
            if(CPsite_search_step < 4){
                CPsite_search_step <- 4
                chr.abun <- proximalAdj(chr.abun, MINSIZE, PolyA_PWM, 
                                        genome, classifier, classifier_cutoff, 
                                        shift_range, search_point_START)
                if(!is.null(tmpfolder)) 
                    save(list=c("chr.abun", "CPsite_search_step"),
                         file=file.path(tmpfolder, seqn))
                if(!silence) 
                    message("chromsome", seqn, "proximal adjust ... done")
            }
        }
        
        chr.abun <- polishCPs(chr.abun)
    }else{
        NULL
    }
}