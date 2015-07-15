prepare4GSEA <- function(eset, groupList, Preranked=TRUE,
                         folder=".", 
                         rnkFilename="InPAS.rnk",
                         chipFilename="InPAS.chip", 
                         dataFilename="dPDUI.txt",
                         PhenFilename="group.cls"){
    if(class(eset)!="UTR3eSet"){
        stop("eset must be an object of UTR3eSet")
    }
    if(Preranked){
        if(length(eset$testRes)<1){
            stop("There is no p.value for ranking. Please try testUsage first")
        }
        testRes <- eset$testRes
        usage <- eset$usage
        if(!identical(usage$transcript, rownames(testRes))){
            stop("The rownames of eset slots are not identical")
        }
        rank <- data.frame(symbol=usage$symbol, p.value=1-testRes[, "P.Value"])
        write.table(rank, file.path(foler, rnkFilename),
                    sep="\t", quote=FALSE, row.names=F, col.names=F)
    }else{
        eset <- as(eset, "ExpressionSet")
        exprs <- exprs(eset)
        samples <- colnames(exprs)
        fD <- fData(eset)
        chip <- fD[, c("transcript", "symbol", "gene")]
        colnames(chip) <- c("Probe Set ID", "Gene Symbol", "Gene Title")
        dir.create(folder, showWarnings=FALSE)
        write.table(chip, file.path(folder, chipFilename), 
                    sep="\t", quote=FALSE, row.names=F)
        exprs <- cbind(NAME=rownames(exprs), exprs)
        write.table(exprs, file.path(folder, dataFilename), 
                    sep="\t", quote=FALSE, row.names=F)
        if((!missing(groupList)) && is.list(groupList)){
            gp <- rep(names(groupList), sapply(groupList, length))
            names(gp) <- unlist(groupList)
            if(all(samples %in% names(gp))){
                gp <- gp[samples]
                levels <- unique(gp)
                out <- c(paste(length(gp), length(levels), 1),
                         paste("#", paste(levels, collapse=" ")),
                         paste(gp, collapse=" "))
                writeLines(out, con=file.path(folder, PhenFilename))
            }
        }
    }
}