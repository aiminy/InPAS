get.regions.coverage <- function(chr, utr3.regions.chr, 
                                 hugeData, coverage, phmm=FALSE){
    view <- utr3.regions.chr[[chr]]
    end <- end(view)
    maxEnd <- max(end)
    if(hugeData){
        .cov <- list()
        if(phmm) all.tx <- list()
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            .ele <- cvg[[chr]]
            if(maxEnd>length(.ele)){
                .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
            }
            if(class(.ele)=="Rle"){
                .cvg <- Views(.ele, start(view), end(view))
                if(phmm) all.tx[[names(coverage)[i]]] <- 
                    viewApply(.cvg, as.integer)
                .cvg <- viewMeans(.cvg, na.rm=TRUE)
            }else{
                ##sum(.ele)==0
                if(sum(.ele)!=0){
                    save.image(file="error.rds")
                    stop("sum of current chromosome is not zero.")
                }
                .cvg <- rep(0, length(view))
            }
            .cov[[names(coverage)[i]]] <- .cvg
            rm(cvg)
        }
    }else{
        .cov <- lapply(coverage, function(.ele){
            .ele <- .ele[[chr]]
            if(maxEnd>length(.ele)){
                .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
            }
            .cvg <- Views(.ele, start(view), end(view))
            .cvg <- viewMeans(.cvg, na.rm=TRUE)
        })
        if(phmm){
            all.tx <- lapply(coverage, function(.ele){
                .ele <- .ele[[chr]]
                if(maxEnd>length(.ele)){
                    .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
                }
                .cvg <- Views(.ele, start(view), end(view))
                viewApply(.cvg, as.integer)
            })
        }
    }
    this.trans <- list()
    for(i in 1:length(.cov[[1]])){
        this.trans[[i]] <- list()
        for(j in names(.cov)){
            this.trans[[i]][[j]] <- .cov[[j]][[i]]
        }
        this.trans[[i]] <- do.call(cbind, this.trans[[i]])
    }
    view$data <- this.trans
    if(phmm) view$data2 <- all.tx[[1]]
    view
}