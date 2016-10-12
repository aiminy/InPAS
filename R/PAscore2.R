PAscore2 <- function(seqname, pos, str, idx, idx.gp,
                     genome, classifier, classifier_cutoff){
    if(length(pos)<1){
        return(NULL)
    }
    coor <- paste(seqname, pos, str, sep="_")
    gr <- 
        GRanges(seqname, 
                IRanges(pos, pos, names=coor), 
                strand=str)
    gr$id <- 1:length(gr)
    coor.id <- !duplicated(coor)
    gr$duplicated <- gr$id[coor.id][match(coor, coor[coor.id])]
    gr.s <- gr[coor.id]
    pred.prob.test <- do.call(rbind, 
                              lapply(split(gr.s, 
                                           rep(1:ceiling(length(gr.s)/1000), 
                                               each=1000)[1:length(gr.s)]), 
                                            function(.gr.s){
        testSet.NaiveBayes <- 
            buildFeatureVector(.gr.s, BSgenomeName = genome, 
                               upstream = classifier@info@upstream,
                               downstream = classifier@info@downstream, 
                               wordSize = classifier@info@wordSize, 
                               alphabet=classifier@info@alphabet,
                               sampleType = "unknown",replaceNAdistance = 30, 
                               method = "NaiveBayes", ZeroBasedIndex = 1, 
                               fetchSeq = TRUE)
        suppressMessages(.pred.prob.test <- 
                             predictTestSet(
                                 testSet.NaiveBayes=testSet.NaiveBayes,
                                 classifier=classifier,
                                 outputFile=NULL, 
                                 assignmentCutoff=classifier_cutoff))
        .pred.prob.test
    }))
    rownames(pred.prob.test) <- NULL
    pred.prob.test <- pred.prob.test[match(names(gr.s), 
                                           pred.prob.test$PeakName), , 
                                     drop=FALSE]
    if(any(duplicated(coor))){
        ## need to recover the order of inputs
        pred.prob.test <- pred.prob.test[match(gr$duplicated, gr.s$id), , 
                                         drop=FALSE]
        pred.prob.test[, "PeakName"] <- names(gr)
    }
    pred.prob.test <- cbind(pred.prob.test[,1:4], idx, idx.gp)
    pred.prob.test <- pred.prob.test[!is.na(pred.prob.test[, "pred.class"]), ]
    pred.prob.test <- pred.prob.test[pred.prob.test[, "pred.class"]==1, ]
    pred.prob.test <- pred.prob.test[order(pred.prob.test[, "idx.gp"],
                                           -pred.prob.test[, "prob True"]), ]
    pred.prob.test <- pred.prob.test[!duplicated(pred.prob.test[, "idx.gp"]), ]
    pred.prob.test
}
