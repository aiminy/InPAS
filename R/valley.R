

## calcultate p value:
## T=abs(mean(sqrt(short_coverage1))/mean(sqrt(long_coverage2)) 
##    - mean(sqrt(long_coverage1))/mean(sqrt(long_coverage2))) 
## sqrt also could be log2, sqrt --> poisson distribution
## log2 --> power law
## p.value = n(Ti<T)/n.total(T)

valley <- function(x, ss, se, n=1, savedID=NA, filterByPval=TRUE){
    pos <- which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == 2) + ss - 1
    if(length(pos)>0 && filterByPval){
        pos.pos <- 
            which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == -2) + ss - 1
        pos.pos <- sort(c(pos, pos.pos, ss, se))
        w <- which(pos.pos %in% pos)
        space=10
        y <- mapply(function(a, b, p) {
            c <- unique(sort(seq(a, b)))
            c <- c[c!=p & c>ss & c<se]
            d <- x[c]
        }, pos.pos[w-1] - space, pos.pos[w-1] + space, pos, SIMPLIFY=FALSE)
        mu <- sapply(y, mean, na.rm=TRUE)
        sigma <- sapply(y, sd, na.rm=TRUE)
        z <- (x[pos]-mu)/sigma
        p <- 2*pnorm(-abs(z))
        pos <- pos[p < 0.001]
    }
    #pos <- unique(c(ss, pos, se))
    pos <- pos[!is.na(x[pos])]
    tobeadd <-ifelse(x[ss] < x[se], ss, se)
    if(length(pos)>0){
        w <- x[tobeadd] < min(x[pos])
        if(length(w)>0 && !is.na(w) && is.logical(w)){
            if(w) pos <- unique(c(pos, tobeadd))
        }
    }else{
        pos <- tobeadd
    }
    if(!is.na(savedID) && savedID>ss && savedID<se) 
        pos <- unique(c(pos, savedID))
    pos <- pos[x[pos]<=quantile(x[ss:se], probs=.5)]
    if(n==-1 || n>length(pos)) n <- length(pos)
    if(length(pos)<1){
        return(pos)
    }else{
        pos <- pos[order(x[pos])][1:n]
    }
}