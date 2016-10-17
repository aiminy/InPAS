inPAS <- function(bedgraphs, genome, utr3, txdb=NA,
                  tags, hugeData=FALSE, ...,
                  
                  gp1, gp2,
                  
                  window_size=100, 
                  search_point_START=50, search_point_END=NA, 
                  cutStart=window_size, cutEnd=0, 
                  coverage_threshold=5, long_coverage_threshold=2, 
                  background=c("same_as_long_coverage_threshold", 
                               "1K", "5K", "10K", "50K"),
                  adjust_distal_polyA_end=TRUE, 
                  PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, 
                  shift_range=window_size, 
                  
                  method=c("limma", "fisher.exact", 
                           "singleSample", "singleGroup"),
                  normalize=c("none", "quantiles", "quantiles.robust",
                              "mean", "median"),
                  design, contrast.matrix, coef=1,
                  
                  P.Value_cutoff=0.05,
                  adj.P.Val_cutoff=0.05, 
                  dPDUI_cutoff=0.3, 
                  PDUI_logFC_cutoff=0.59,
                  
                  BPPARAM=NULL){
    background <- match.arg(background)
    method <- match.arg(method)
    normalize <- match.arg(normalize)
    
    if(!missing(gp1) | !missing(gp2)){
        if(!all(c(gp1,gp2) %in% tags)) stop("gp1 and gp2 must be in tags")
        if(missing(gp2)) gp2 <- NULL
        if(missing(gp1)) gp1 <- NULL
        groupList <- list(gp1=gp1, gp2=gp2)
    }else{
        groupList <- NULL
    }
    
    if(missing(gp2)) gp2 <- NULL
    if(missing(gp1)) gp1 <- NULL
    if(length(gp2)>0) gp2 <- gp2[!is.na(gp2)]
    if(length(gp1)>0) gp1 <- gp1[!is.na(gp1)]
    if(missing(design)) design <- NULL
    if(missing(contrast.matrix)) contrast.matrix <- NULL
    
    if(!missing(txdb)){
        if(background!="same_as_long_coverage_threshold") {
            stop("txdb is missing when you want local background")
        }else{
            if(class(txdb)!="TxDb")
                stop("txdb must be an object of TxDb")
        }
    }
    
    if(missing(genome) || missing(utr3))
        stop("genome and utr3 are required.")
    if(class(genome)!="BSgenome")
        stop("genome must be an object of BSgenome.")
    if(class(utr3)!="GRanges" | 
           !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))){
        stop("utr3 must be output of function of utr3Annotation")
    }
    
    if(length(gp2)<1 | length(gp1)<1){
        samples <- unlist(groupList)
        samples <- samples[!is.na(samples)]
        if(length(samples)>1){
            if(method!="singleGroup")
                stop("method should be singleGroup")
        }else{
            if(length(samples)==1){
                if(method!="singleSample")
                    stop("method should be singleSample")
            }
        }
    }
    
    ##step1 coverage
    coverage <- 
        coverageFromBedGraph(bedgraphs, tags, genome, hugeData=hugeData, ...)
    
    ##step2 predict CPsites
    CPs <- CPsites(coverage=coverage, groupList=groupList, 
                   genome=genome, utr3=utr3, 
                   window_size=window_size, 
                   search_point_START=search_point_START, 
                   search_point_END=search_point_END, 
                   cutStart=cutStart, cutEnd=cutEnd, 
                   adjust_distal_polyA_end=adjust_distal_polyA_end, 
                   background=background,
                   txdb=txdb,
                   coverage_threshold=coverage_threshold, 
                   long_coverage_threshold=long_coverage_threshold,
                   PolyA_PWM=PolyA_PWM, 
                   classifier=classifier, 
                   classifier_cutoff=classifier_cutoff, 
                   shift_range=shift_range, BPPARAM=BPPARAM)
    
    ##step3 calculate usage
    res <- testUsage(CPsites=CPs, 
                     coverage=coverage, 
                     genome=genome,
                     utr3=utr3, 
                     method=method,
                     normalize=normalize,
                     design=design, 
                     contrast.matrix=contrast.matrix,
                     coef=coef,
                     gp1=gp1, gp2=gp2)
    
    if(hugeData){
        removeTmpfile(coverage)
    }
    
    filterRes(res, gp1=gp1, gp2=gp2,
              background_coverage_threshold=long_coverage_threshold,
              P.Value_cutoff=P.Value_cutoff,
              adj.P.Val_cutoff=adj.P.Val_cutoff, 
              dPDUI_cutoff=dPDUI_cutoff, 
              PDUI_logFC_cutoff=PDUI_logFC_cutoff)
}

removeTmpfile <- function(coverage){
    lapply(coverage, unlink)
}