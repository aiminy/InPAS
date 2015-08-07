CPsites <- function(coverage, groupList=NULL, genome, utr3, window_size=100, 
                    search_point_START=50, search_point_END=NA, 
                    cutStart=window_size, cutEnd=0, adjust_distal_polyA_end=TRUE, 
                    coverage_threshold=5, long_coverage_threshold=2, 
                    background=c("same_as_long_coverage_threshold", 
                                 "1K", "5K", "10K", "50K"),
                    txdb=NA,
                    gcCompensation=NA, mappabilityCompensation=NA, 
                    FFT=FALSE, fft.sm.power=20, 
                    PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, 
                    shift_range=window_size, BPPARAM=NULL){
    if(!is.na(PolyA_PWM)[1]){
        if(class(PolyA_PWM)!="matrix") stop("PolyA_PWM must be matrix")
        if(any(rownames(PolyA_PWM)!=c("A", "C", "G", "T"))) 
            stop("rownames of PolyA_PWM must be c(\"A\", \"C\", \"G\", \"T\")")
    }
    if(missing(coverage) || missing(genome) || missing(utr3))
        stop("coverage, genome and utr3 are required.")
    if(class(genome)!="BSgenome")
        stop("genome must be an object of BSgenome.")
    if(class(utr3)!="GRanges" | 
           !all(utr3$id %in% c("utr3", "next.exon.gap", "CDS"))){
        stop("utr3 must be output of function of utr3Annotation")
    }
    if(seqlevelsStyle(utr3)!=seqlevelsStyle(genome)){
        stop("the seqlevelsStyle of utr3 must be same as genome")
    }
    background <- match.arg(background)
    introns <- GRanges()
    if(background!="same_as_long_coverage_threshold"){
        if(class(txdb)!="TxDb") stop("txdb is missing when you want local background")
        if(!identical(unname(genome(txdb))[1], unname(genome(utr3))[1])) 
            stop("genome of txdb should be same as genome of utr3")
        introns <- intronsByTranscript(txdb, use.names=TRUE)
        introns <- unlist(introns)
        exons <- exons(txdb)
        exons <- disjoin(exons, ignore.strand=TRUE)
        introns <- disjoin(introns, ignore.strand=TRUE)
        intron_exons <- disjoin(c(exons, introns))
        ol.exons <- findOverlaps(exons, intron_exons)
        introns <- intron_exons[-unique(subjectHits(ol.exons))]
    }
    utr3 <- utr3[utr3$id!="CDS"]
    MINSIZE <- 10
    hugeData <- class(coverage[[1]])=="character"
    if(hugeData && !is.null(names(groupList))[1]){
        coverage <- mergeCoverage(coverage, groupList, genome, BPPARAM=BPPARAM)
        hugeData <- FALSE
        depth.weight <- depthWeight(coverage, hugeData)
        totalCov <- totalCoverage(coverage, genome, hugeData)
    }else{
        depth.weight <- depthWeight(coverage, hugeData, groupList)
        totalCov <- totalCoverage(coverage, genome, hugeData, groupList)
    }
    utr3TotalCov <- UTR3TotalCoverage(utr3, totalCov, 
                                      gcCompensation, mappabilityCompensation, 
                                      FFT=FFT, fft.sm.power=fft.sm.power)
    z2s <- zScoreThrethold(background, introns, totalCov, utr3)
    if(!is.null(BPPARAM)){
        shorten_UTR_estimation <- 
            bplapply(utr3TotalCov, CPsite_estimation,
                     BPPARAM=BPPARAM, utr3=utr3,
                     MINSIZE=MINSIZE, 
                     window_size=window_size, 
                     search_point_START=search_point_START,
                     search_point_END=search_point_END, 
                     cutStart=cutStart, cutEnd=cutEnd, 
                     adjust_distal_polyA_end=adjust_distal_polyA_end, 
                     background=background,
                     z2s=z2s,
                     coverage_threshold=coverage_threshold, 
                     long_coverage_threshold=long_coverage_threshold, 
                     PolyA_PWM=PolyA_PWM, 
                     classifier=classifier, 
                     classifier_cutoff=classifier_cutoff, 
                     shift_range=shift_range, 
                     depth.weight=depth.weight, 
                     genome=genome)
    }else{
        shorten_UTR_estimation <- 
            lapply(utr3TotalCov, CPsite_estimation, 
                   utr3=utr3, MINSIZE=MINSIZE, 
                   window_size=window_size, 
                   search_point_START=search_point_START, 
                   search_point_END=search_point_END, 
                   cutStart=cutStart, cutEnd=cutEnd, 
                   adjust_distal_polyA_end=adjust_distal_polyA_end,  
                   background=background,
                   z2s=z2s,
                   coverage_threshold=coverage_threshold, 
                   long_coverage_threshold=long_coverage_threshold, 
                   PolyA_PWM=PolyA_PWM, 
                   classifier=classifier, classifier_cutoff=classifier_cutoff, 
                   shift_range=shift_range, 
                   depth.weight=depth.weight, 
                   genome=genome)
    }
    shorten_UTR_estimation <- 
        do.call(rbind, 
                unname(shorten_UTR_estimation[!sapply(shorten_UTR_estimation, 
                                                      is.null)]))
    utr3.shorten.UTR <- utr3[utr3$id=="utr3"]
    utr3.shorten.UTR$id <- NULL
    utr3.shorten.UTR <- 
        utr3.shorten.UTR[utr3.shorten.UTR$transcript 
                         %in% rownames(shorten_UTR_estimation)]
    shorten_UTR_estimation <- 
        shorten_UTR_estimation[utr3.shorten.UTR$transcript, , drop=FALSE]
    utr3.shorten.UTR$fit_value <- unlist(shorten_UTR_estimation[,"fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <- 
        unlist(shorten_UTR_estimation[,"Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <- 
        unlist(shorten_UTR_estimation[,"Predicted_Distal_APA"])
    utr3.shorten.UTR$type <- unlist(shorten_UTR_estimation[,"type"])
    utr3.shorten.UTR$utr3start <- unlist(shorten_UTR_estimation[,"utr3start"])
    utr3.shorten.UTR$utr3end <- unlist(shorten_UTR_estimation[,"utr3end"])
    utr3.shorten.UTR <- 
        utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
}