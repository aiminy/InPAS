coverageRate <- function(coverage, txdb, genome, 
                         cutoff_readsNum=1,
                         cutoff_expdGene_cvgRate=0.1,
                         cutoff_expdGene_sampleRate=0.5,
                         which=NULL, ...){
    stopifnot(class(txdb)=="TxDb")
    hugeData <- class(coverage[[1]])=="character"
    stopifnot(class(genome)=="BSgenome")
    seqnames <- trimSeqnames(genome)
    if(length(which)>0){
        stopifnot(class(which)=="GRanges")
        message("strand information will be ignored.")
        seqnames <- seqnames[seqnames %in% unique(as.character(seqnames(which)))]
    }
    seqLen <- seqLen(genome)
    seqRle <- sapply(seqLen, function(.ele) Rle(0, .ele), simplify = FALSE)
    
    exon <- exons(txdb, columns=c("exon_id", "gene_id"))
    exon$gene_id <- sapply(exon$gene_id, `[`, 1)
    exon_gene_id_map <- exon$gene_id
    names(exon_gene_id_map) <- as.character(exon$exon_id)
    UTR3 <- threeUTRsByTranscript(txdb)
    UTR3 <- unlist(UTR3, use.names=FALSE)
    stopifnot(all(UTR3$exon_id %in% exon$exon_id))
    UTR3$gene_id <- exon_gene_id_map[as.character(UTR3$exon_id)]
    mcols(UTR3) <- mcols(UTR3)[, c("exon_id", "gene_id")]
    
    exon <- exon[!is.na(exon$gene_id)]
    UTR3 <- UTR3[!is.na(UTR3$gene_id)] 
    
    
    if(length(which)>0){
        ol <- findOverlaps(exon, which, ignore.strand=TRUE)
        exon <- exon[sort(unique(queryHits(ol)))]
        ol <- findOverlaps(UTR3, which, ignore.strand=TRUE)
        UTR3 <- UTR3[sort(unique(queryHits(ol)))]
    }
    
    ## reduce exon for each gene
    reduce_by_gene <- function(gr){
        stopifnot(class(gr)=="GRanges")
        stopifnot(!is.null(gr$gene_id))
        seqn <- cbind(as.character(seqnames(gr)), gr$gene_id)
        seqn <- unique(seqn)
        seqn_map <- seqn[, 1]
        names(seqn_map) <- seqn[, 2]
        gr.fake <- GRanges(seqnames=gr$gene_id, 
                             ranges=ranges(gr),
                             strand=strand(gr))
        gr.fake <- reduce(gr.fake)
        GRanges(seqn_map[as.character(seqnames(gr.fake))], 
                ranges=ranges(gr.fake),
                strand=strand(gr.fake),
                gene_id=as.character(seqnames(gr.fake)))
    }
    
    exon <- reduce_by_gene(exon)
    UTR3 <- reduce_by_gene(UTR3)
    
    feature.cvg <- function(cov, feature, seqnames, seqRle){
        if(hugeData){
            cvg <- NULL
            load(cov)
            cov <- cvg
        }
        cov <- cov[seqnames]
        names(cov) <- seqnames
        idx <- sapply(cov, is.null)
        cov[idx] <- seqRle[idx]
        cvg.base <- mapply(function(cov.seq, feature.seq){
            vw <- Views(cov.seq, start=start(feature.seq), end=end(feature.seq))
            viewApply(vw, function(.ele){
                sum(unlist(.ele>cutoff_readsNum))
            })
        }, cov[seqnames], 
        feature[seqnames], 
        SIMPLIFY = FALSE)
        stopifnot(identical(sapply(cvg.base, length), sapply(feature, length)))
        unlist(cvg.base, use.names = FALSE)
    }
    
    exon.s <- split(exon, as.character(seqnames(exon)))
    seqnames <- seqnames[seqnames %in% names(exon.s)]
    if(length(seqnames)<1){
        stop("No chromosome is selected.")
    }
    exon.s <- exon.s[seqnames]
    seqRle <- seqRle[seqnames]
    exon.unlist <- 
        if(class(exon.s)!="GRangesList") GRangesList(exon.s) else exon.s
    exon.unlist <- unlist(exon.s)
    
    exon.cvg <- lapply(coverage, feature.cvg, feature=exon.s, 
                       seqnames=seqnames, seqRle=seqRle)
    names(exon.cvg) <- names(coverage)
    exon.cvg <- do.call(cbind, exon.cvg)
    exon.cvg.wid <- cbind(exon.cvg, gene.width=width(exon.unlist))
    gene.cvg <- rowsum(exon.cvg.wid, exon.unlist$gene_id, reorder=FALSE)
    gene.expd.rate <- gene.cvg[, -ncol(gene.cvg)] / gene.cvg[, ncol(gene.cvg)]
    gene.expd.rate <- gene.expd.rate > cutoff_expdGene_cvgRate
    gene.expd.overall <- rowSums(gene.expd.rate) > 
        (cutoff_expdGene_sampleRate * ncol(gene.expd.rate))
    gene.expd.overall.rate <- mean(gene.expd.overall)
    gene.cvg.rate <- colSums(gene.cvg[, -ncol(gene.cvg)]) / sum(gene.cvg[, ncol(gene.cvg)])
    if(gene.expd.overall.rate>0){
        gene.cvg.rate.subsetBy.expd.gene <- 
            colSums(gene.cvg[gene.expd.overall, -ncol(gene.cvg)]) / 
            sum(gene.cvg[gene.expd.overall, ncol(gene.cvg)])
    }else{
        gene.cvg.rate.subsetBy.expd.gene <- rep(0, length(coverage))
        names(gene.cvg.rate.subsetBy.expd.gene) <- names(coverage)
    }
    
    
    UTR3.s <- split(UTR3, as.character(seqnames(UTR3)))
    seqnames <- seqnames[seqnames %in% names(UTR3.s)]
    if(length(seqnames)<1){
        stop("No chromosome is selected.")
    }
    UTR3.s <- UTR3.s[seqnames]
    seqRle <- seqRle[seqnames]
    UTR3.unlist <- if(class(UTR3.s)!="GRangesList") GRangesList(UTR3.s) else UTR3.s
    UTR3.unlist <- unlist(UTR3.s)
    
    UTR3.cvg <- lapply(coverage, feature.cvg, feature=UTR3.s, 
                       seqnames=seqnames, seqRle=seqRle)
    names(UTR3.cvg) <- names(coverage)
    UTR3.cvg <- do.call(cbind, UTR3.cvg)
    UTR3.cvg.rate <- colSums(UTR3.cvg) / sum(width(UTR3.unlist))
    idx <- UTR3.unlist$gene_id %in% names(gene.expd.overall[gene.expd.overall])
    if(sum(idx)>0){
        UTR3.cvg.rate.subsetBy.expd.gene <- 
            colSums(UTR3.cvg[idx, , drop=FALSE]) / sum(width(UTR3.unlist[idx]))
    }else{
        UTR3.cvg.rate.subsetBy.expd.gene <- rep(0, length(coverage))
        names(UTR3.cvg.rate.subsetBy.expd.gene) <- names(coverage)
    }
    
    data.frame(gene.coverage.rate=gene.cvg.rate,
               expressed.gene.coverage.rate=gene.cvg.rate.subsetBy.expd.gene,
               UTR3.coverage.rate=UTR3.cvg.rate,
               UTR3.expressed.gene.subset.coverage.rate=UTR3.cvg.rate.subsetBy.expd.gene)
}