test_CPsites_utr3Usage<-function(){
    utr3 <- GRanges("chr1", 
                    IRanges(c(100, 600), c(599, 2000),
                            names=c("transcript1|gene1|utr3", 
                                    "transcript1|gene1|next.exon.gap")), 
                    strand="+", feature="unknown", 
                    id=c("utr3", "next.exon.gap"),
                    exon="transcript1", transcript="transcript1|",
                    gene="1", symbol="gene1")
    testGR <- GRanges("chr1", IRanges(c(5, 400), c(399, 1000)),
                      strand="*", score=c(40, 10))
    filename <- tempfile()
    export(testGR, filename, format="BEDGraph")
    genome <- BSgenome.Mmusculus.UCSC.mm10
    
    coverage <- coverageFromBedGraph(filename, tags="test", 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, gp1="test", gp2=NULL, genome=genome, 
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5,
                  )
    checkEquals(CP[1]$Predicted_Proximal_APA, 399)
    checkEquals(CP[1]$Predicted_Distal_APA, 999)
    checkEquals(CP[1]$type, "novalDistal")
    res <- utr3UsageEstimation(CP, coverage, gp1="test", gp2=NULL, 
                               short_coverage_threshold=10, 
                               long_coverage_threshold=5)
    
    testGR2 <- GRanges("chr1", IRanges(5, 1000), strand="*", score=20)
    filename2 <- tempfile()
    export(testGR2, filename2, format="BEDGraph")
    coverage <- coverageFromBedGraph(c(filename, filename2), 
                                     tags=c("test1", "test2"), 
                                     genome=genome, 
                                     hugeData=FALSE)
    CP <- CPsites(coverage=coverage, gp1="test1", gp2="test2", genome=genome, 
                  utr3=utr3, coverage_threshold=5, long_coverage_threshold=5,
    )
    checkEquals(CP[1]$Predicted_Proximal_APA, 399)
    checkEquals(CP[1]$Predicted_Distal_APA, 999)
    checkEquals(CP[1]$type, "novalDistal")
    
    res <- utr3UsageEstimation(CP, coverage, gp1="test1", gp2="test2", 
                               short_coverage_threshold=10, 
                               long_coverage_threshold=5)
    
    checkEqualsNumeric(res$PDUI.gp1, 0.25, tolerance=1.0e-2)
    checkEqualsNumeric(res$PDUI.gp2, 1, tolerance=1.0e-2)
    checkEquals(as.logical(res$filterPass), TRUE)
}