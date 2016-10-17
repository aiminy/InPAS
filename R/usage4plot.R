usage4plot <- function(gr, coverage, proximalSites, genome, 
                       groupList=NULL){ 
    gcCompensation=NA
    mappabilityCompensation=NA 
    FFT=FALSE
    fft.sm.power=20
    if(class(gr)!="GRanges") stop("gr must be an object of GRanges")
    if(is.null(names(gr))) 
        names(gr) <- 
            paste("X", formatC(1:length(gr), 
                               width=nchar(length(gr)), flag="0"), sep="")
    if(is.null(names(proximalSites))){
        names(proximalSites) <- names(gr)
    }else{
        if(!all(names(gr) %in% names(proximalSites))) 
            stop("not all GRanges has proximalSites")
    }
    
    hugeData <- class(coverage[[1]])=="character"
    depth.weight <- depthWeight(coverage, hugeData)
    totalCov <- totalCoverage(coverage, genome, hugeData, groupList)
    ## get coverage for each region
    utr3TotalCov <- 
        UTR3TotalCoverage(gr, totalCov, 
                          gcCompensation=gcCompensation, 
                          mappabilityCompensation=mappabilityCompensation, 
                          FFT=FFT, fft.sm.power=fft.sm.power)
    ## get least square error
    grs <- split(gr, as.character(seqnames(gr)))
    seqnames <- names(grs)
    datInfo <- lapply(utr3TotalCov[seqnames], function(.cov){
        .gr <- gr[names(.cov)]
        data <- mapply(function(.ele, .str){
            if(hugeData){
                if(.str) .ele <- rev(.ele)
                se <- length(.ele)-1
                os <- optimalSegmentation(.ele, 
                                          search_point_START=1, 
                                          search_point_END=se)
                cov_diff <- os$cov_diff
                cvg <- .ele[1:se]
            }else{
                if(.str) .ele <- .ele[nrow(.ele):1,, drop=FALSE]
                se <- nrow(.ele)-1
                os <- apply(.ele, 2, optimalSegmentation, 
                            search_point_START=1, search_point_END=se)
                cov_diff <- sapply(os, "[[", "cov_diff")
                cov_diff <- rowMeans(cov_diff)
                cvg <- .ele[1:se, , drop=FALSE]
            }
            cbind(cov_diff, cvg)
        }, .cov, as.character(strand(.gr))=="-", SIMPLIFY=FALSE)
        .gr$dat <- data
        .gr$offset <- ifelse(strand(.gr)=="+", 
                             proximalSites[names(.cov)] - start(.gr),
                             end(.gr) - proximalSites[names(.cov)] + 1)
        .gr
    })
    datInfo <- unlist(GRangesList(datInfo))
}