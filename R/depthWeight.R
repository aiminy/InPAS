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
            groups <- rep(names(groupList), sapply(groupList, length))
            names(groups) <- unlist(groupList)
            groups <- groups[match(names(depth), names(groups))]
            if(any(is.na(groups))){
                stop("all the tags in groupList must has correct names as it in coverage")
            }
            depth <- split(depth, groups)
            depth <- sapply(depth, mean)
            n <- names(depth)
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