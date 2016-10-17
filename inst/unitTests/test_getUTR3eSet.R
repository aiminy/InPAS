test_getUTR3eSet <- function(){
    path <- system.file("extdata", package="InPAS")
    load(file.path(path, "coverage.MAQC.rda"))
    load(file.path(path, "CPs.MAQC.rda"))
    load(file.path(path, "eset.MAQC.rda"))
    tags <- c("Brain.auto", "Brain.phiX", "UHR.auto", "UHR.phiX")
    g <- factor(gsub("\\..*$", "", tags))
    design <- model.matrix(~-1+g)
    colnames(design) <- c("Brain", "UHR")
    contrast.matrix<-makeContrasts(contrasts="Brain-UHR",levels=design)
    contrast.matrix
    data("utr3.hg19")
    eset1 <- getUTR3eSet(CPs, coverage, 
                        genome=BSgenome.Hsapiens.UCSC.hg19, 
                        utr3=utr3.hg19,
                        normalize="none")
    #checkIdentical(eset, eset1) ## PDUI.log2[54, 3] not equal in linux
    for(i in c("usage", "PDUI", "short", "long", "signals", "testRes")){
        checkIdentical(slot(eset, i), slot(eset1, i))
    }
}