test_coverageFromBedGraph<-function(){
    start <- seq(1, 9, by=2)
    end <- seq(2, 10, by=2)
    testGR <- GRanges("chr1", IRanges(start, end), strand="+", score=1:5)
    filename <- tempfile()
    export(testGR, filename, format="BEDGraph")
    cov <- coverageFromBedGraph(filename, tags="test", 
                                BSgenome.Mmusculus.UCSC.mm10, hugeData=FALSE)
    checkEquals(rep(2, 5), runLength(cov$test$chr1)[1:5])
    checkEquals(1:5, runValue(cov$test$chr1)[1:5])
}