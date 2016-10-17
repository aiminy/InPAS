##remove utr3---___---utr3, need to improve.
removeUTR3__UTR3 <- function(x){
    ## smooth by window_size=10
    ws <- 10
    len <- length(x)
    if(len>100){
        y <- rowsum(x, group=rep(1:ceiling(len/ws), each=ws)[1:len], 
                    reorder=FALSE)
        y <- y[-length(y)]
        id <- valley(y, 1, length(y), 1, filterByPval=FALSE)
        if(length(id)>0){
            if(id<3) id <- 3
            pval <- try(t.test(y[1:(id-2)], y[-(1:(id+2))])$p.value, 
                        silent=TRUE)
            if(class(pval)=="numeric" && length(pval)==1 && !is.na(pval)){
                if(pval < 0.001){
                    x <- x[1:((id-1)*ws)]
                }
            }
        }
    }
    x
}