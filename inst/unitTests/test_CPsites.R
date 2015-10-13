checkCPs <- function(CP, proximal, distal, type){
    if(!missing(proximal)) 
        checkEquals(CP[1]$Predicted_Proximal_APA, proximal)
    if(!missing(distal))
        checkEquals(CP[1]$Predicted_Distal_APA, distal)
    if(!missing(type))
        checkEquals(CP[1]$type, type)
}

generateGRs <- function(starts, ends, scores){
    GRanges("chr1", IRanges(start=starts, end=ends), 
            strand="*", score=scores)
}

generateUTR3 <- function(normal=TRUE){
    if(normal){
        GRanges("chr1", 
                IRanges(c(100, 600), c(599, 2000),
                        names=c("transcript1|gene1|utr3", 
                                "transcript1|gene1|next.exon.gap")), 
                strand="+", feature="unknown", 
                id=c("utr3", "next.exon.gap"),
                exon="transcript1", transcript="transcript1|",
                gene="1", symbol="gene1")
    }else{
        data(utr3.mm10)
        utr3.mm10[seqnames(utr3.mm10)=="chr1" & end(utr3.mm10)<10000000]
    }
}

NOtest_searchDistalCPs <- function(){
    # drop from end
    utr3 <- generateUTR3()
    testGR <- generateGRs(c(5, 400, 1001), c(399, 1000, 1500), c(40, 10, 1))
    
    filename <- tempfile()
    export(testGR, filename, format="BEDGraph")
    genome <- BSgenome.Mmusculus.UCSC.mm10
    
    coverage <- coverageFromBedGraph(filename, tags="test", 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, genome=genome, 
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5)
    checkCPs(CP, 399, 1000, "novel distal")
    
    # local background
    utr3 <- generateUTR3(FALSE)
    #     signal, background=2, signal/background=2.5 for uc007afk.2_5|Rgs20, 
    #     GRanges object with 3 ranges and 1 metadata column:
    #         seqnames             ranges strand |     score
    #            <Rle>          <IRanges>  <Rle> | <numeric>
    #     [1]     chr1 [4909200, 4909475]      * |         5
    #     [2]     chr1 [4909476, 4910473]      * |        10
    #     [3]     chr1 [4910474, 4910662]      * |        11
    
    background <- rep(rnbinom(140000, 10e9, mu=2), each=50)
    signal <- rep(c(0, 5, 10, 11, 0), c(1909199, 276, 998, 189, 5089338)) +
        background
    signal <- signal[1899200:1920662]
    signal <- rle(signal)
    signal.cumsum <- cumsum(signal$lengths)
    testGR <- generateGRs(4899200+c(0, signal.cumsum[-length(signal.cumsum)]), 
                          4899199+signal.cumsum,
                          signal$values)
    export(testGR, filename, format="BEDGraph")
    genome <- BSgenome.Mmusculus.UCSC.mm10
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    coverage <- coverageFromBedGraph(filename, tags="test", 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, genome=genome, txdb=txdb,
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5,
                  background="10K")
    checkCPs(CP, proximal=4909476, type="novel distal")
    
    ## distal adjust
#     data(classifier) ## could not load cleanUpdTSeq from the test file. need to figure out the reason.
#     CP <- CPsites(coverage=coverage, genome=genome, txdb=txdb,
#                   utr3=utr3, coverage_threshold=5, long_coverage_threshold=5,
#                   background="10K", classifier=classifier, 
#                   classifier_cutoff=.8)
    
    ## remove utr3---___---utr3,
    
}

test_CPsites_utr3Usage<-function(){
    utr3 <- generateUTR3()
    testGR <- generateGRs(c(5, 400), c(399, 1000), score=c(40, 10))
    filename <- tempfile()
    export(testGR, filename, format="BEDGraph")
    genome <- BSgenome.Mmusculus.UCSC.mm10
    
    coverage <- coverageFromBedGraph(filename, tags="test", 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, genome=genome, 
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5)
    checkCPs(CP, 399, 1000, "novel distal")
    
    res <- utr3UsageEstimation(CP, coverage, genome=genome, 
                               utr3=utr3, gp1="test", gp2=NULL, 
                               short_coverage_threshold=10, 
                               long_coverage_threshold=5)
    
    testGR2 <- GRanges("chr1", IRanges(5, 1000), strand="*", score=20)
    filename2 <- tempfile()
    export(testGR2, filename2, format="BEDGraph")
    coverage <- coverageFromBedGraph(c(filename, filename2), 
                                     tags=c("test1", "test2"), 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, genome=genome, 
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5)
    checkCPs(CP, 399, 1000, "novel distal")
    
    res <- utr3UsageEstimation(CP, coverage, genome=genome, 
                               utr3=utr3, gp1="test1", gp2="test2", 
                               short_coverage_threshold=10, 
                               long_coverage_threshold=5)
    
    checkEqualsNumeric(res$PDUI.gp1, 0.25, tolerance=1.0e-2)
    checkEqualsNumeric(res$PDUI.gp2, 1, tolerance=1.0e-2)
    checkEquals(as.logical(res$filterPass), TRUE)
}

test_getCov_seqlevelsStyle <- function(){
    utr3 <- generateUTR3()
    testGR <- generateGRs(c(5, 400, 1001), c(399, 1000, 1500), c(40, 10, 1))
    seqlevelsStyle(testGR) <- "NCBI"
    
    filename <- tempfile()
    export(testGR, filename, format="BEDGraph")
    genome <- BSgenome.Mmusculus.UCSC.mm10
    
    coverage <- coverageFromBedGraph(filename, tags="test", 
                                     genome=genome, 
                                     hugeData=FALSE)
}