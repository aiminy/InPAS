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