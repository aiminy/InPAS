PAscore <- function(seqname, pos, str, idx, PWM, genome, ups=50, dws=50){
    pos <- pos[!is.na(pos)]
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