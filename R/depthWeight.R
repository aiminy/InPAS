depthWeight <- function(coverage, hugeData, groupList=NULL){
    n <- names(coverage)
    if(hugeData){
        depth <- numeric(length(coverage))
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            d <- sapply(cvg, function(.cvg) {
                sum(as.double(runValue(.cvg)) * runLength(.cvg))
            })
            depth[i] <- sum(d)
            rm(cvg)
        }
        if(!is.null(names(groupList))[1]){
            names(depth) <- names(coverage)
            groups <- do.call(rbind, groupList)
            depth <- split(depth[groups[,1]], rownames(groups))
            depth <- lapply(depth, sum)
            n <- rownames(groups)
        }
        
    }else{
        depth <- sapply(coverage, function(cvg){
            d <- sapply(cvg, function(.cvg) {
                sum(as.double(runValue(.cvg)) * runLength(.cvg))
            })
            sum(d)
        })
        
    }
    
    depth.weight <- depth / mean(depth)
    names(depth.weight) <- n
    
    depth.weight
}