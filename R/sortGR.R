sortGR <- function(.ele){
    str <- as.character(strand(.ele))[1] == "+"
    if(str){ ## sort the element by position
        .ele <- .ele[order(start(.ele))]
    }else{
        .ele <- .ele[order(-end(.ele))]
    }
    .ele
}