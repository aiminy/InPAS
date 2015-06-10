totalCoverage <- function(coverage, genome, hugeData, groupList=NULL){
    if(!hugeData) return(coverage)
    ## calculate total coverage
    seqnames <- trimSeqnames(genome)
    seqLen <- seqLen(genome)
    
    if(!is.null(names(groupList))[1]){
        cov <- vector("list", length=length(groupList))
        names(cov) <- names(groupList)
        for(i in 1:length(cov)){
            for(s in seqnames){
                cov[[i]][[s]] <- Rle(0, seqLen[s])
            }
        }
        
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            gp <- names(coverage)[i]
            gp <- which(sapply(groupList, function(.ele) any(gp %in% .ele)))[1]
            if(!is.null(gp)){
                for(s in seqnames){
                    if(s %in% names(cvg)){
                        cov[[gp]][[s]] <- cov[[gp]][[s]] + cvg[[s]]
                    }
                }
            }
            rm(cvg)
        }
        idx <- rep(FALSE, length(seqnames))
        names(idx) <- seqnames
        for(i in 1:length(cov)){
            for(s in seqnames){
                if(nrun(cov[[i]][[s]])!=1 || runValue(cov[[i]][[s]][1])!=0) 
                    idx[s] <- TRUE
            }
        }
        for(i in 1:length(cov)){
            cov[[i]] <- cov[[i]][idx]
        }
    }else{
        cov <- vector("list", length=length(coverage))
        names(cov) <- names(coverage)
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            gp <- names(coverage)[i]
            cov[[gp]] <- list()
            for(s in seqnames){
                if(s %in% names(cvg)){
                    cov[[gp]][[s]] <- cvg[[s]]
                }
            }
            rm(cvg)
        }
    }
    cov
}