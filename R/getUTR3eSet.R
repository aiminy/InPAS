getUTR3eSet <- function(CPsites, coverage, genome, utr3, 
                        normalize=c("none", "quantiles", "quantiles.robust",
                                 "mean", "median"),
                        ...,
                        BPPARAM=NULL, singleSample=FALSE){
    if(missing(coverage) || missing(CPsites))
        stop("CPsites and coverage is required.")
    if(missing(utr3) || missing(genome)){
        stop("utr3 and genome is required.")
    }
    if(class(genome)!="BSgenome")
        stop("genome must be an object of BSgenome.")
    if(class(utr3)!="GRanges" | 
           !all(utr3$id %in% c("utr3", "next.exon.gap", "CDS"))){
        stop("utr3 must be output of function of utr3Annotation")
    }
    normalize <- match.arg(normalize)
    hugeData <- class(coverage[[1]])=="character"
    if((!singleSample) && length(coverage)==1){
        message("Single sample mode is on")
        singleSample <- TRUE
    }
    if(singleSample){
        if(length(coverage)>1) message("Only first sample will be used.")
        UTRusage <- UTR3usage(CPsites, coverage, hugeData, BPPARAM, phmm=TRUE)
    }else{
        UTRusage <- UTR3usage(CPsites, coverage, hugeData, BPPARAM)
    }
    
    UTRusage <- split(UTRusage, UTRusage$transcript)
    UTRusage <- UTRusage[sapply(UTRusage, length)==2]
    UTRusage <- unlist(UTRusage, use.names = FALSE)
    UTRusage.total <- UTRusage[UTRusage$source=="short"]
    UTRusage.long <- UTRusage[UTRusage$source=="long"]
    UTRusage.total <- UTRusage.total[!duplicated(UTRusage.total$transcript)]
    UTRusage.long <- UTRusage.long[match(UTRusage.total$transcript, 
                                         UTRusage.long$transcript)]
    PDUItable <- UTRusage.total
    PDUItable$data <- NULL
    start(PDUItable)[as.character(strand(PDUItable))=="-"] <- 
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable))=="-"]
    end(PDUItable)[as.character(strand(PDUItable))=="+"] <- 
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable))=="+"]
    
    
    UTRusage.total.data <- do.call(rbind, UTRusage.total$data)
    UTRusage.long.data <- do.call(rbind, UTRusage.long$data)
    UTRusage.short.data <- UTRusage.total.data - UTRusage.long.data
    lt0 <- apply(UTRusage.short.data, 1, function(.ele) any(.ele<0))
    if(any(lt0)){
        CDS <- utr3[utr3$id=="CDS"]
        CDS <- CDS[CDS$transcript %in% unique(PDUItable$transcript[lt0])]
        CDSusage <- lastCDSusage(CDS, coverage, hugeData, BPPARAM)
        CDSusage.data <- do.call(rbind, CDSusage$data)
        rownames(CDSusage.data) <- CDSusage$transcript
        idx <- match(PDUItable$transcript[lt0], CDSusage$transcript)
        
        UTRusage.short.data[lt0[!is.na(idx)]] <- 
            CDSusage.data[idx[!is.na(idx)],] - UTRusage.long.data[lt0[!is.na(idx)]]
        UTRusage.short.data[UTRusage.short.data<0] <- 0
    }
    ##normalization?
    normalize.foo <- function(exprs, FUN){
        avgs = FUN(exprs)$Estimates
        scaling.factors = avgs /avgs[1]
        scaling.factors = matrix(rep(scaling.factors, nrow(exprs)), 
                                 ncol=length(scaling.factors), byrow=TRUE)
        exprs = exprs /scaling.factors
    }
    normalize.mean <- function(exprs){
        normalize.foo(exprs, colSummarizeAvg)
    }
    normalize.median <- function(exprs){
        normalize.foo(exprs, colSummarizeMedian)
    }
    if(normalize!="none"){
        UTRusage.long.short.data <- rbind(UTRusage.long.data, 
                                         UTRusage.short.data)
        UTRusage.long.short.data <- 
            switch(normalize,
                   quantile=normalize.quantiles(UTRusage.long.short.data),
                   quantile.robust=normalize.quantiles.robust(UTRusage.long.short.data, ...),
                   mean=normalize.mean(UTRusage.long.short.data),
                   median=normalize.median(UTRusage.long.short.data),
                   UTRusage.long.short.data
            )## suppose the output should be same order
        UTRusage.long.data <- 
            UTRusage.long.short.data[1:nrow(UTRusage.long.data), , drop=FALSE]
        UTRusage.short.data <- 
            UTRusage.long.short.data[-(1:nrow(UTRusage.long.data)), , drop=FALSE]
    }
    
    UTRusage.PDUI <- UTRusage.long.data/(UTRusage.long.data + UTRusage.short.data)
    UTRusage.PDUI.log2 <- log2(UTRusage.PDUI+.Machine$double.xmin)
    rownames(UTRusage.PDUI) <- 
        rownames(UTRusage.PDUI.log2) <- 
        rownames(UTRusage.long.data) <- 
        rownames(UTRusage.short.data) <- PDUItable$transcript
    PDUItable$source <- NULL
    if(singleSample){
        signals.short <- UTRusage.total$data2
        signals.long <- UTRusage.long$data2
        cut50 <- function(x, y){
            z <- c(x, y)
            if(length(z)>50) {
                xat <- floor(50*length(x)/length(z))
                z <- tapply(z, cut(1:length(z), 50), mean)
            }else{
                xat <- length(x)
            }
            c(xat, z)
        }
        signals <- mapply(cut50, signals.short, signals.long, SIMPLIFY=FALSE)
        names(signals) <- UTRusage.total$transcript
        new("UTR3eSet", usage=PDUItable, PDUI=UTRusage.PDUI,
             PDUI.log2=UTRusage.PDUI.log2, 
             short=UTRusage.short.data,
             long=UTRusage.long.data,
             signals=signals)
    }else{
        new("UTR3eSet", usage=PDUItable, PDUI=UTRusage.PDUI,
             PDUI.log2=UTRusage.PDUI.log2, 
             short=UTRusage.short.data,
             long=UTRusage.long.data)
    }
}