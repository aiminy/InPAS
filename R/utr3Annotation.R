utr3Annotation <- function(txdb, orgDbSYMBOL, MAX_EXONS_GAP=10000){
    if(class(txdb)!="TxDb") stop("txdb must be an object of TxDb")
    if(!grepl("^org.*egSYMBOL$", orgDbSYMBOL)) 
        stop("orgDbSYMBOL must be a string 
             start with org and end with egSYMBOL.")
    exons <- ranges(GeneRegionTrack(txdb))
    exons$density <- NULL
    exons.old <- exons ## need to keep all the exons, for genome fragments
    exons <- exons[exons$feature!="ncRNA"]
    utr3.all <- exons[exons$feature=="utr3"]
    utr3 <- utr3.all[!is.na(utr3.all$gene)]
    tryCatch(env<-get(orgDbSYMBOL),error = function(e){
        stop(paste("Annotation database ",
                   orgDbSYMBOL,
                   " does not exist!\n\t",
                   "Please try to load annotation database by library(",
                   gsub("SYMBOL$", "", orgDbSYMBOL), ".db)", sep=""),
             call.=FALSE)
    })
    symbol <- AnnotationDbi::mget(utr3$gene, env, ifnotfound=NA)
    symbol <- sapply(symbol, paste, collapse=";")
    utr3$symbol <- symbol
    
    ##get last 3utr for each tx
    ## keep the last 3utr for each transcripts
    utr3.plus <- utr3[strand(utr3)=="+"]
    utr3.minus <- utr3[strand(utr3)=="-"]
    utr3.plus <- utr3.plus[order(-start(utr3.plus))]
    utr3.minus <- utr3.minus[order(end(utr3.minus))]
    ###get last 3utr
    utr3.plus <- utr3.plus[!duplicated(utr3.plus$transcript)]
    utr3.minus <- utr3.minus[!duplicated(utr3.minus$transcript)]
    utr3.last <- c(utr3.plus, utr3.minus)
    exons.not.utr3.last <- exons.old[!exons.old$exon %in% utr3.last$exon]
    
    getRange <- function(gr, f){
        seq <- as.character(seqnames(gr))
        st <- start(gr)
        en <- end(gr)
        strand <- as.character(strand(gr))
        f <- as.character(mcols(gr)[,f])
        seq <- split(seq, f)
        st <- split(st, f)
        en <- split(en, f)
        strand <- split(strand, f)
        st <- sapply(st, min)
        en <- sapply(en, max)
        seq <- sapply(seq, `[`, 1)
        strand <- sapply(strand, `[`, 1)
        if((!identical(names(st), names(en))) || 
               (!identical(names(st), names(seq))) || 
               (!identical(names(st), names(strand)))) 
            stop("unexpect happend at getRange!")
        GRanges(seq, IRanges(st, en, names=names(st)), strand)
    }
    
    ## extract utr3 sharing same start positions from same gene 
    if(any(is.na(utr3.last$symbol))) 
        stop("unexpect happend at check gene symbol")
    start.utr3.last <- 
        paste(as.character(seqnames(utr3.last)), start(utr3.last))
    end.utr3.last <- paste(as.character(seqnames(utr3.last)), end(utr3.last))
    start.dup <- 
        unique(start.utr3.last[duplicated(start.utr3.last) & 
                                   as.character(strand(utr3.last))=="+"])
    end.dup <- 
        unique(end.utr3.last[duplicated(end.utr3.last) & 
                                 as.character(strand(utr3.last))=="-"])
    utr3.last$feature <- "unknown"
    utr3.last.dup <- 
        utr3.last[(start.utr3.last %in% start.dup) | 
                      (end.utr3.last %in% end.dup)]
    if(length(utr3.last.dup)>0){
        utr3.last.nd <- 
            utr3.last[!((start.utr3.last %in% start.dup) | 
                            (end.utr3.last %in% end.dup))]
        dup.group.plus <- paste(as.character(seqnames(utr3.last.dup)), 
                                start(utr3.last.dup))
        dup.group.minus <- paste(as.character(seqnames(utr3.last.dup)), 
                                 end(utr3.last.dup))
        dup.group <- ifelse(as.character(strand(utr3.last.dup))=="+", 
                            dup.group.plus, 
                            dup.group.minus)
        utr3.last.dup <- split(utr3.last.dup, dup.group)
        utr3.last.dup <- lapply(utr3.last.dup, function(.ele){
            ## assum all the symbol should be same
            if(length(unique(.ele$symbol))!=1) return(.ele[-(1:length(.ele))])
            str <- as.character(strand(.ele))[1]
            if(str=="+"){
                .ele <- .ele[order(end(.ele))]
            }else{
                .ele <- .ele[order(-start(.ele))]
            }
            wid <- width(.ele)
            dif <- diff(wid)
            .ele <- .ele[c(wid[1]>1, dif>1)]
            len <- length(.ele)
            if(len<=1) return(.ele)
            if(str=="+"){
                .ele$feature <- paste("proximalCP", 
                                      paste(end(.ele)[-len], 
                                            collapse="_"), sep="_")
            }else{
                .ele$feature <- paste("proximalCP", 
                                      paste(start(.ele)[-len], 
                                            collapse="_"), sep="_")
            }
            .ele[len]
        })
        ## attract longest last 3UTR
        utr3.last.dup <- utr3.last.dup[sapply(utr3.last.dup, length)>0]
        utr3.last.dup <- unlist(GRangesList(utr3.last.dup))
        names(utr3.last.dup) <- NULL
        utr3.last <- c(utr3.last.nd, utr3.last.dup)
    }
    utr3.last.block <- utr3.last
    
    ## mask 5UTR-CDS and any other region from last 3UTR
    non.utr3 <- getRange(exons.not.utr3.last, "transcript")
    
    removeMask <- function(gr, mask){
        mask.reduce <- reduce(mask)
        ol.mask <- findOverlaps(mask.reduce, gr, 
                                ignore.strand=TRUE, minoverlap=1)
        if(length(ol.mask)>0){
            gr.ol <- gr[subjectHits(ol.mask)]
            gr.nol <- gr[-unique(subjectHits(ol.mask))]
            mask.ol <- mask.reduce[queryHits(ol.mask)]
            gr.ol$rel <- "ol"
            gr.ol$rel[end(mask.ol)>=start(gr.ol) & 
                          start(mask.ol)<=start(gr.ol)] <- "left"
            gr.ol$rel[end(mask.ol)>=end(gr.ol) & 
                          start(mask.ol)<=end(gr.ol)] <- "right"
            gr.ol$rel[end(mask.ol)>=end(gr.ol) & 
                          start(mask.ol)<=start(gr.ol)] <- "cover"
            end(gr.ol)[gr.ol$rel=="right"] <- 
                start(mask.ol)[gr.ol$rel=="right"]
            start(gr.ol)[gr.ol$rel=="left"] <- end(mask.ol)[gr.ol$rel=="left"]
            start(gr.ol)[gr.ol$rel=="ol" & 
                             as.character(strand(gr.ol))=="+"] <- 
                end(mask.ol)[gr.ol$rel=="ol" & 
                                 as.character(strand(gr.ol))=="+"] + 1
            end(gr.ol)[gr.ol$rel=="ol" & as.character(strand(gr.ol))=="-"] <- 
                start(mask.ol)[gr.ol$rel=="ol" & 
                                   as.character(strand(gr.ol))=="-"] - 1
            gr.ol <- gr.ol[gr.ol$rel!="cover"]
            gr.ol$rel <- NULL
            ##leave last segment for each utr3
            exons.dup <- gr.ol[gr.ol$exon %in% 
                                   gr.ol$exon[duplicated(gr.ol$exon)]]
            if(length(exons.dup)>0){
                exons.nd <- 
                    gr.ol[!gr.ol$exon %in% gr.ol$exon[duplicated(gr.ol$exon)]]
                exons.dup <- split(exons.dup, exons.dup$exon)
                exons.dup <- lapply(exons.dup, function(.ele){
                    ## get overlaps parts
                    .disj <- disjoin(.ele)
                    .cnt <- countOverlaps(.disj, .ele)
                    .disj <- .disj[order(-.cnt)]
                    .disj <- .disj[1]
                    mcols(.disj) <- mcols(.ele[1])
                    .disj
                })
                exons.dup <- unlist(GRangesList(exons.dup))
                names(exons.dup) <- NULL
                if(length(exons.nd)>0) gr.ol <- c(exons.dup, exons.nd) 
                else gr.ol <- exons.dup
            }
            
            gr <- c(gr.ol, gr.nol)
        }else{
            gr
        }
    }
    utr3.last <- removeMask(utr3.last, non.utr3)

    ## cut off overlaps of reverse 3UTR
    utr3.block <- threeUTRsByTranscript(txdb)
    utr3.block <- unlist(utr3.block)
    utr3.block$tx_id <- names(utr3.block)
    utr3.block <- getRange(utr3.block, "tx_id")
    utr3.reverse.strand <- utr3.block
    strand(utr3.reverse.strand) <- 
        ifelse(as.character(strand(utr3.reverse.strand))=="+", "-", "+")
    utr3.reverse.strand <- reduce(utr3.reverse.strand)
    ol.utr3.rev <- findOverlaps(utr3.last, utr3.reverse.strand, 
                                ignore.strand=FALSE, minoverlap=1)
    if(length(ol.utr3.rev)>0){
        utr3.last.ol <- utr3.last[queryHits(ol.utr3.rev)]
        utr3.last.nol <- utr3.last[-unique(queryHits(ol.utr3.rev))]
        utr3.rev.ol <- utr3.reverse.strand[subjectHits(ol.utr3.rev)]
        strand <- as.character(strand(utr3.last.ol))=="+"
        idx <- (start(utr3.last.ol) < start(utr3.rev.ol) & strand) | 
            (end(utr3.last.ol) > end(utr3.rev.ol) & !strand)
        utr3.last.ol <- utr3.last.ol[idx]
        utr3.rev.ol <- utr3.rev.ol[idx]
        strand <- strand[idx]
        end(utr3.last.ol)[strand] <- start(utr3.rev.ol)[strand] - 1
        start(utr3.last.ol)[!strand] <- end(utr3.rev.ol)[!strand] + 1
        utr3.last <- c(utr3.last.nol, utr3.last.ol[width(utr3.last.ol)>1])
    }
        
    ## keep the last 3utr for each transcripts
    utr3.last.plus <- utr3.last[strand(utr3.last)=="+"]
    utr3.last.minus <- utr3.last[strand(utr3.last)=="-"]
    utr3.last.plus <- utr3.last.plus[order(-start(utr3.last.plus))]
    utr3.last.minus <- utr3.last.minus[order(end(utr3.last.minus))]
    ###get last 3utr
    utr3.last.plus <- utr3.last.plus[!duplicated(utr3.last.plus$transcript)]
    utr3.last.minus <- utr3.last.minus[!duplicated(utr3.last.minus$transcript)]
    utr3.last <- c(utr3.last.plus, utr3.last.minus)
    
    ## remove all the utr3 belong to single exons
    single.exon.transcripts <- split(exons.old$feature, exons.old$transcript)
    single.exon.transcripts <- 
        sapply(single.exon.transcripts, function(.ele) sum(.ele=="CDS"))
    single.exon.transcripts <- 
        names(single.exon.transcripts)[single.exon.transcripts<1]
    utr3.last <- utr3.last[!utr3.last$transcript %in% single.exon.transcripts]
    
    ## remove all the exactly same 3utr for same gene.
    utr3.last <- unique(utr3.last) 
    utr3.last <- utr3.last[width(utr3.last)>1]
    
    
    ## cut left side for overlaps of left same strand 3UTR.
    ## find last sigment for disjoin of block
    utr3.last.disjoin <- 
        disjoin(c(disjoin(utr3.last), disjoin(utr3.last.block)))
    utr3.last.disjoin <- utr3.last.disjoin[width(utr3.last.disjoin)>1]
    ol <- findOverlaps(utr3.last, utr3.last.disjoin, 
                       ignore.strand=FALSE, minoverlap=1L)
    ol.query <- utr3.last[queryHits(ol)]
    ol.subject <- utr3.last.disjoin[subjectHits(ol)]
    mcols(ol.subject) <- mcols(ol.query)
    ## get last of utr3 for each transcript
    utr3.last.plus <- ol.subject[strand(ol.subject)=="+"]
    utr3.last.minus <- ol.subject[strand(ol.subject)=="-"]
    utr3.last.plus <- utr3.last.plus[order(-start(utr3.last.plus))]
    utr3.last.minus <- utr3.last.minus[order(end(utr3.last.minus))]
    ###get last 3utr
    utr3.last.plus <- utr3.last.plus[!duplicated(utr3.last.plus$transcript)]
    utr3.last.minus <- utr3.last.minus[!duplicated(utr3.last.minus$transcript)]
    utr3.last <- c(utr3.last.plus, utr3.last.minus)
    ## handle same start
    utr3.last.sameStart <- utr3.last[grepl("proximalCP", utr3.last$feature)]
    utr3.last.unknown <- utr3.last[utr3.last$feature=="unknown"]
    utr3.last.sameStart.CP <- 
        strsplit(gsub("proximalCP_", "", utr3.last.sameStart$feature), "_")
    utr3.last.sameStart.Start <- start(utr3.last.sameStart)
    utr3.last.sameStart.End <- end(utr3.last.sameStart)
    utr3.last.sameStart$feature <- mapply(function(cp, st, en){
        cp <- as.integer(cp)
        cp <- cp[cp>st & cp<en]
        if(length(cp)==0){
            return("unknown")
        }else{
            return(paste("proximalCP", paste(cp, collapse="_"), sep="_"))
        }
    }, utr3.last.sameStart.CP, 
    utr3.last.sameStart.Start, 
    utr3.last.sameStart.End, 
    SIMPLIFY=TRUE)
    utr3.clean <- c(utr3.last.sameStart, utr3.last.unknown)
    utr3.clean <- 
        utr3.clean[order(as.character(seqnames(utr3.clean)), 
                         start(utr3.clean))]
    
    ## get extented utr3 region
    ## get follow exons
    exons.rd <- exons.old
    strand(exons.rd) <- "*"
    exons.rd <- reduce(exons.rd)
    gaps <- gaps(exons.rd)
    utr3.clean.ext1 <- utr3.clean
    str <- as.character(strand(utr3.clean.ext1))=="+"
    utr3.clean.ext1[str] <- 
        shift(utr3.clean.ext1, shift=1, use.names=FALSE)[str]
    utr3.clean.ext1[!str] <- 
        shift(utr3.clean.ext1, shift=-1, use.names=FALSE)[!str]
    ol <- findOverlaps(utr3.clean.ext1, gaps, ignore.strand=TRUE)
    ol.utr3.clean <- utr3.clean[queryHits(ol)]
    next.exons.gap <- gaps[subjectHits(ol)]
    strand(next.exons.gap) <- strand(ol.utr3.clean)
    mcols(next.exons.gap) <- mcols(ol.utr3.clean)
    wid <- width(next.exons.gap)>MAX_EXONS_GAP
    width(next.exons.gap)[wid & as.character(strand(next.exons.gap))=="+"] <- 
        MAX_EXONS_GAP
    start(next.exons.gap)[wid & as.character(strand(next.exons.gap))=="-"] <- 
        end(next.exons.gap)[
            wid & as.character(strand(next.exons.gap))=="-"] - MAX_EXONS_GAP
    ## lable the next.exons as special label
    next.exons.gap$id <- "next.exon.gap"
    utr3.clean$id <- "utr3"
    
    utr3.fixed <- c(utr3.clean, next.exons.gap)
    names(utr3.fixed) <- 
        paste(utr3.fixed$exon, utr3.fixed$symbol, utr3.fixed$id, sep="|")
    utr3.fixed
}