mergeCoverage <- function(coverage, groupList, genome, BPPARAM=NULL){
    seqnames <- trimSeqnames(genome)
    seqLen <- seqLen(genome)
    seqRle <- list()
    for(i in 1:length(seqnames)){
        seqRle[[seqnames[i]]] <- Rle(0, seqLen[i])
    }
    
    gpl <- rep(names(groupList), sapply(groupList, length))
    names(gpl) <- unlist(groupList)
    coverageL <- split(coverage[names(gpl)], gpl)
    
    cov <- list()
    
    for(j in 1:length(coverageL)){
        covFiles <- coverageL[[j]]
        x <- 1:length(covFiles)
        y <- split(x, ceiling(x/10))
        
        cov[[j]] <- seqRle
        
        for(i in 1:length(y)){
            if(!is.null(BPPARAM)){
                cv <- bptry(bplapply(covFiles[y[[i]]], function(.ele){
                    cvg <- NULL
                    load(.ele)
                    cvg <- cvg[seqnames]
                    names(cvg) <- seqnames
                    idx <- sapply(cvg, is.null)
                    cvg[idx] <- seqRle[idx]
                    cvg
                }, BPPARAM=BPPARAM))
                while(!all(bpok(cv))){
                    cv <- bptry(bplapply(covFiles[y[[i]]], function(.ele){
                        cvg <- NULL
                        load(.ele)
                        cvg <- cvg[seqnames]
                        names(cvg) <- seqnames
                        idx <- sapply(cvg, is.null)
                        cvg[idx] <- seqRle[idx]
                        cvg
                    }, BPREDO=cv, BPPARAM=BPPARAM))
                }
                if(length(cv)<10){
                    for(m in (length(cv)+1):10){
                        cv[[m]] <- seqRle
                    }
                }
                cov[[j]] <- bptry(bplapply(seqnames, function(s){
                    cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
                        cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
                        cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
                        cv[[10]][[s]] + cov[[j]][[s]]
                }, BPPARAM=BPPARAM))
                while(!all(bpok(cov[[j]]))){
                    cov[[j]] <- bptry(bplapply(seqnames, function(s){
                        cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
                            cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
                            cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
                            cv[[10]][[s]] + cov[[j]][[s]]
                    }, BPREDO=cov[[j]], BPPARAM=BPPARAM))
                }
            }else{
                cv <- lapply(covFiles[y[[i]]], function(.ele){
                    cvg <- NULL
                    load(.ele)
                    cvg <- cvg[seqnames]
                    names(cvg) <- seqnames
                    idx <- sapply(cvg, is.null)
                    cvg[idx] <- seqRle[idx]
                    cvg
                })
                if(length(cv)<10){
                    for(m in (length(cv)+1):10){
                        cv[[m]] <- seqRle
                    }
                }
                cov[[j]] <- lapply(seqnames, function(s){
                    cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
                        cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
                        cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
                        cv[[10]][[s]] + cov[[j]][[s]]
                })
            }
            names(cov[[j]]) <- seqnames
        }
    }
    names(cov) <- names(coverageL)
    cov
}