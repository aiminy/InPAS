getUTR3region <- function(.grs){
    .ra <- c(ranges(.grs), 
             IRanges(min(.grs$Predicted_Proximal_APA[1], 
                         .grs$Predicted_Distal_APA[1]), 
                     max(.grs$Predicted_Proximal_APA[1], 
                         .grs$Predicted_Distal_APA[1]), 
                     names="CPsite"))
    .ra <- disjoin(.ra)
    .ra <- .ra[width(.ra)>1]
    if(as.character(strand(.grs))[1]=="+"){
        short <- .ra[start(.ra)<.grs$Predicted_Proximal_APA[1], ]
        long <- .ra[start(.ra)>=.grs$Predicted_Proximal_APA[1], ]
    }else{
        long <- .ra[start(.ra)<.grs$Predicted_Proximal_APA[1], ]
        short <- .ra[start(.ra)>=.grs$Predicted_Proximal_APA[1], ]
    }
    long <- reduce(long)
    short <- reduce(short)
    if(length(long)==0){
        gr <- GRanges(as.character(seqnames(.grs))[1], short, 
                      strand=as.character(strand(.grs))[1], 
                      source=rep("short", length(short)))
    }else{
        if(length(short)==0){
            gr <- GRanges(as.character(seqnames(.grs))[1], long, 
                          strand=as.character(strand(.grs))[1], 
                          source=rep("long", length(long)))
        }else{
            gr <- c(GRanges(as.character(seqnames(.grs))[1], long, 
                            strand=as.character(strand(.grs))[1], 
                            source=rep("long", length(long))),
                    GRanges(as.character(seqnames(.grs))[1], short, 
                            strand=as.character(strand(.grs))[1], 
                            source=rep("short", length(short))))
        }
    }
#     gr$transcript <- .grs$transcript[1]
#     gr$gene <- .grs$gene[1]
#     gr$symbol <- .grs$symbol[1]
#     gr$fit_value <- .grs$fit_value[1]
#     gr$Predicted_Proximal_APA <- .grs$Predicted_Proximal_APA[1]
#     gr$Predicted_Distal_APA <- .grs$Predicted_Distal_APA[1]
#     gr$type <- .grs$type[1]
#     gr$utr3start <- .grs$utr3start[1]
#     gr$utr3end <- .grs$utr3end[1]
    for(mc in colnames(mcols(.grs))){
        if(mc!="exon"){
            mcols(gr)[, mc] <- mcols(.grs)[1, mc]
        }
    }
    gr
}