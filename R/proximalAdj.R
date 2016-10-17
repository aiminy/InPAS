proximalAdj <- function(CPs, MINSIZE, PolyA_PWM, 
                        genome, classifier, classifier_cutoff, 
                        shift_range, search_point_START, step=1){
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
                                              shift_range, search_point_START,
                                              step)
    }
    CPs$Predicted_Proximal_APA[flag] <- idx.list[flag]
    
    CPs
}