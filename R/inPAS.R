inPAS <- function(bedgraphs, tags, genome, utr3, gp1, gp2, txdb, 
                  window_size=100, short_coverage_threshold=NA, 
                  long_coverage_threshold=NA, 
                  adjusted.P_val.cutoff=0.05, 
                  dPDUI_cutoff=0.3, 
                  PDUI_logFC_cutoff=0.59, 
                  search_point_START=50, search_point_END=NA, 
                  cutStart=100, cutEnd=0, 
                  search_distal_polyA_end=FALSE, 
                  coverage_threshold=5, 
                  gcCompensation=NA, 
                  mappabilityCompensation=NA, 
                  FFT=FALSE, fft.sm.power=20, 
                  hugeData=FALSE, 
                  PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, 
                  shift_range=0, BPPARAM=NULL){
    if(!all(c(gp1,gp2) %in% tags)) stop("gp1 and gp2 must be in tags")
    ##step1 coverage
    coverage <- 
        coverageFromBedGraph(bedgraphs, tags, genome, hugeData=hugeData)
    
    if(is.na(long_coverage_threshold)||is.na(short_coverage_threshold)){
        covT <- covThreshold(coverage, genome, txdb, utr3, 
                             hugeData=hugeData, gp1=gp1, gp2=gp2)
        if(is.na(long_coverage_threshold)){
            long_coverage_threshold <- covT[1]
        }
        if(is.na(short_coverage_threshold)){
            short_coverage_threshold <- covT[2]
        }
    }
    ##step2 predict CPsites
    CPsites <- CPsites(coverage=coverage, gp1=gp1, gp2=gp2, 
                       genome=genome, utr3=utr3, 
                       window_size=window_size, 
                       search_point_START=search_point_START, 
                       search_point_END=search_point_END, 
                       cutStart=cutStart, cutEnd=cutEnd, 
                       search_distal_polyA_end=search_distal_polyA_end, 
                       coverage_threshold=coverage_threshold, 
                       long_coverage_threshold=long_coverage_threshold, 
                       gcCompensation=gcCompensation, 
                       mappabilityCompensation=mappabilityCompensation, 
                       FFT=FFT, fft.sm.power=fft.sm.power, 
                       PolyA_PWM=PolyA_PWM, 
                       classifier=classifier, 
                       classifier_cutoff=classifier_cutoff, 
                       shift_range=shift_range, BPPARAM=BPPARAM)
    ##step3 calculate usage
    res <- 
        utr3UsageEstimation(CPsites, coverage, gp1, gp2, 
                            short_coverage_threshold=short_coverage_threshold, 
                            long_coverage_threshold=long_coverage_threshold, 
                            adjusted.P_val.cutoff=adjusted.P_val.cutoff, 
                            dPDUI_cutoff=dPDUI_cutoff, 
                            PDUI_logFC_cutoff=PDUI_logFC_cutoff, 
                            BPPARAM=BPPARAM)
    if(hugeData){
        removeTmpfile(coverage)
    }
    res
}

removeTmpfile <- function(coverage){
    lapply(coverage, unlink)
}