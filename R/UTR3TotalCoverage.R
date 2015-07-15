compensation <- function(view, comp, start, end){
    mapply(function(.ele, .s, .e){
        .ele * comp[.s:.e]
    }, view, start, end, SIMPLIFY=FALSE)
}


##smooth methods
fft.smooth <- function(sn, p){#
    if(length(sn)<=p) return(sn)
    sn.fft <- fft(sn)#
    sn.fft[p:length(sn.fft)] <- 0+0i#
    sn.ifft = fft(sn.fft, inverse = TRUE)/length(sn.fft)#
    Re(sn.ifft)#
}

UTR3TotalCoverage <- function(utr3, totalCov, 
                              gcCompensation=NA, 
                              mappabilityCompensation=NA, 
                              FFT=FALSE, 
                              fft.sm.power=20){
    utr3.s <- split(utr3, as.character(seqnames(utr3)))
    seqnames <- unlist(sapply(totalCov, names))
    seqnames <- table(seqnames)
    seqnames <- names(seqnames)[seqnames==length(totalCov)]
    seqnames <- sort(intersect(names(utr3.s), seqnames))
    cov <- list()
    for(seq in seqnames){
        cov[[seq]] <- list()
        for(n in names(totalCov)){
            cov[[seq]][[n]] <- totalCov[[n]][[seq]]
        }
    }
    mapply(function(.utr, .cvg, gcComp, mapComp){ ##memory consume
        strand <- as.character(strand(.utr))
        start <- start(.utr)
        end <- end(.utr)
        maxEnd <- max(end)
        for(i in 1:length(.cvg)){
            if(maxEnd>length(.cvg[[i]])){
                .cvg[[i]] <- 
                    append(.cvg[[i]], rep(0, maxEnd - length(.cvg[[i]])+1))
            }
        }
        
        view <- lapply(.cvg, Views, start=start, end=end, names=names(.utr))
        view <- sapply(view, function(.view){
            viewApply(.view, as.integer)
        })
        if(class(view[1,])=="list"){
            view <- apply(view, 1, function(.ele) do.call(cbind, .ele))
        }else{
            view <- list(view)
            names(view) <- names(.utr)
        }
        if(length(view)>0) view <- view[sapply(view, mean)>0]
        if(!is.na(gcComp[1]) && length(gcComp)==length(.cvg[[1]])){
            view <- compensation(view, gcComp, start, end)
        }
        if(!is.na(mapComp[1]) && length(mapComp)==length(.cvg[[1]])){
            view <- compensation(view, mapComp, start, end)
        }
        if(FFT){
            lapply(view, function(.ele) 
                apply(.ele, 2, fft.smooth, p=fft.sm.power))
        }else{
            view
        }
    }, utr3.s[seqnames], cov[seqnames], 
    gcCompensation[seqnames], 
    mappabilityCompensation[seqnames], 
    SIMPLIFY=FALSE)
}