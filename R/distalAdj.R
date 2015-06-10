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