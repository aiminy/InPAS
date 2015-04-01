getUTR3region <- function(.grs){
    .ra <- c(ranges(.grs), 
             IRanges(min(.grs$Predicted_Proximal_APA[1], 
                         .grs$Predicted_Distal_APA[1]), 
                     max(.grs$Predicted_Proximal_APA[1], 
                         .grs$Predicted_Distal_APA[1]), 
                     names="CPsite"))
    .ra <- disjoin(.ra)
    .ra <- .ra[width(.ra)>1]
    if(as.character(strand(.grs))[1]=="+"){
        short <- .ra[start(.ra)<.grs$Predicted_Proximal_APA[1], ]
        long <- .ra[start(.ra)>=.grs$Predicted_Proximal_APA[1], ]
    }else{
        long <- .ra[start(.ra)<.grs$Predicted_Proximal_APA[1], ]
        short <- .ra[start(.ra)>=.grs$Predicted_Proximal_APA[1], ]
    }
    long <- reduce(long)
    short <- reduce(short)
    if(length(long)==0){
        gr <- GRanges(as.character(seqnames(.grs))[1], short, 
                      strand=as.character(strand(.grs))[1], 
                      source=rep("short", length(short)))
    }else{
        if(length(short)==0){
            gr <- GRanges(as.character(seqnames(.grs))[1], long, 
                          strand=as.character(strand(.grs))[1], 
                          source=rep("long", length(long)))
        }else{
            gr <- c(GRanges(as.character(seqnames(.grs))[1], long, 
                            strand=as.character(strand(.grs))[1], 
                            source=rep("long", length(long))),
                    GRanges(as.character(seqnames(.grs))[1], short, 
                            strand=as.character(strand(.grs))[1], 
                            source=rep("short", length(short))))
        }
    }
    gr$transcript <- .grs$transcript[1]
    gr$gene <- .grs$gene[1]
    gr$symbol <- .grs$symbol[1]
    gr$fit_value <- .grs$fit_value[1]
    gr$Predicted_Proximal_APA <- .grs$Predicted_Proximal_APA[1]
    gr$Predicted_Distal_APA <- .grs$Predicted_Distal_APA[1]
    gr$type <- .grs$type[1]
    gr$utr3start <- .grs$utr3start[1]
    gr$utr3end <- .grs$utr3end[1]
    gr
}
get.utr3.regions.coverage <- function(chr, utr3.regions.chr, 
                                      hugeData, coverage, phmm=FALSE){
    view <- utr3.regions.chr[[chr]]
    end <- end(view)
    maxEnd <- max(end)
    if(hugeData){
        .cov <- list()
        if(phmm) all.tx <- list()
        for(i in 1:length(coverage)){
            cvg <- NULL
            load(coverage[[i]])
            .ele <- cvg[[chr]]
            if(maxEnd>length(.ele)){
                .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
            }
            .cvg <- Views(.ele, start(view), end(view))
            if(phmm) all.tx[[names(coverage)[i]]] <- 
                viewApply(.cvg, as.integer)
            .cvg <- viewMeans(.cvg, na.rm=TRUE)
            .cov[[names(coverage)[i]]] <- .cvg
            rm(cvg)
        }
    }else{
        .cov <- lapply(coverage, function(.ele){
            .ele <- .ele[[chr]]
            if(maxEnd>length(.ele)){
                .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
            }
            .cvg <- Views(.ele, start(view), end(view))
            .cvg <- viewMeans(.cvg, na.rm=TRUE)
        })
        if(phmm){
            all.tx <- lapply(coverage, function(.ele){
                .ele <- .ele[[chr]]
                if(maxEnd>length(.ele)){
                    .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
                }
                .cvg <- Views(.ele, start(view), end(view))
                viewApply(.cvg, as.integer)
            })
        }
    }
    this.trans <- list()
    for(i in 1:length(.cov[[1]])){
        this.trans[[i]] <- list()
        for(j in names(.cov)){
            this.trans[[i]][[j]] <- .cov[[j]][[i]]
        }
        this.trans[[i]] <- do.call(cbind, this.trans[[i]])
    }
    view$data <- this.trans
    if(phmm) view$data2 <- all.tx[[1]]
    view
}
UTR3usage <- function(CPsites, coverage, hugeData, BPPARAM=NULL, phmm=FALSE){
    utr3.selected.shorten.UTR <- 
        split(CPsites, as.character(seqnames(CPsites)))
    seqnames <- names(utr3.selected.shorten.UTR)
    
    ##covert coverage from list of Rle to list of matrix
    ##step1 prepare the short and long utr3 regions
    utr3.trans.shorten.UTR <- split(CPsites, CPsites$transcript)
    if(!is.null(BPPARAM)){
        utr3.regions <- bplapply(utr3.trans.shorten.UTR, 
                                 getUTR3region, BPPARAM=BPPARAM)
    }else{
        utr3.regions <- lapply(utr3.trans.shorten.UTR, getUTR3region)
    }
    utr3.regions <- unlist(GRangesList(utr3.regions))
    ##step2, calculate coverage and merge into a matrix for long and short
    utr3.regions.chr <- split(utr3.regions, seqnames(utr3.regions))
    utr3.regions.chr <- utr3.regions.chr[seqnames]
    if(!is.null(BPPARAM)){
        utr3.regions.chr <- bplapply(seqnames, get.utr3.regions.coverage, 
                                     BPPARAM=BPPARAM, 
                                     utr3.regions.chr=utr3.regions.chr, 
                                     hugeData=hugeData, 
                                     coverage=coverage, 
                                     phmm=phmm)
    }else{
        utr3.regions.chr <- lapply(seqnames, get.utr3.regions.coverage, 
                                   utr3.regions.chr=utr3.regions.chr, 
                                   hugeData=hugeData, 
                                   coverage=coverage, 
                                   phmm=phmm)
    }
    cov.ranges <- unlist(GRangesList(utr3.regions.chr))
}

utr3UsageEstimation <- function(CPsites, coverage, gp1, gp2=NULL, 
                                short_coverage_threshold=10, 
                                long_coverage_threshold=2, 
                                adjusted.P_val.cutoff=0.05, 
                                dPDUI_cutoff=0.3, 
                                PDUI_logFC_cutoff=0.59, 
                                BPPARAM=NULL){
    if(!all(c(gp1,gp2) %in% names(coverage))) 
        stop("gp1 and gp2 must be in names of coverage")
    hugeData <- class(coverage[[1]])=="character"
    if(length(c(gp1, gp2))==1){
        coverage <- coverage[c(gp1, gp2)]
        UTRusage <- UTR3usage(CPsites, coverage, hugeData, BPPARAM, phmm=TRUE)
    }else{
        depth.weight <- depthWeight(coverage, hugeData)
        UTRusage <- UTR3usage(CPsites, coverage, hugeData, BPPARAM)
    }
    ##step3, calculate mean for each group
    UTRusage <- split(UTRusage, UTRusage$transcript)
    UTRusage <- UTRusage[sapply(UTRusage, length)>=2]
    UTRusage <- unlist(UTRusage, use.names = FALSE)
    UTRusage.short <- UTRusage[UTRusage$source=="short"]
    UTRusage.long <- UTRusage[UTRusage$source=="long"]
    UTRusage.short <- UTRusage.short[!duplicated(UTRusage.short$transcript)]
    UTRusage.long <- UTRusage.long[match(UTRusage.short$transcript, 
                                         UTRusage.long$transcript)]
    PDUItable <- UTRusage.short
    PDUItable$data <- NULL
    start(PDUItable)[as.character(strand(PDUItable))=="-"] <- 
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable))=="-"]
    end(PDUItable)[as.character(strand(PDUItable))=="+"] <- 
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable))=="+"]
    
    UTRusage.short.data <- do.call(rbind, UTRusage.short$data)
    UTRusage.long.data <- do.call(rbind, UTRusage.long$data)
    PDUItable$total.gp1 <- UTRusage.short.data[, gp1]
    PDUItable$long.gp1 <- UTRusage.long.data[, gp1]
    PDUItable$short.gp1 <- PDUItable$total.gp1 - PDUItable$long.gp1
    PDUItable$short.gp1[PDUItable$short.gp1<0] <- 0
    if(length(c(gp1, gp2))==1){
        PDUItable$data2 <- NULL
        PDUItable$test_status <- PDUItable$total.gp1>long_coverage_threshold & 
            PDUItable$long.gp1>long_coverage_threshold
        ##test by phmm ##phmm not ready, simple as t.test
        pval <- mapply(function(a, b){
            if(length(a)>10) a <- tapply(a, cut(1:length(a), 10), sum)
            if(length(b)>10) b <- tapply(b, cut(1:length(b), 10), sum)
            t.test(a, b)$p.value
        }, UTRusage.short$data2, UTRusage.long$data2)
        PDUItable$dPDUI <- PDUItable$PDUI.gp1 <- 
            PDUItable$long.gp1/(PDUItable$long.gp1 + PDUItable$short.gp1)
        PDUItable$logFC <- TRUE
    }else{
        if(is.character(gp2)){
            PDUItable$total.gp2 <- UTRusage.short.data[, gp2]
            PDUItable$long.gp2 <- UTRusage.long.data[, gp2]
            PDUItable$short.gp2 <- PDUItable$total.gp2 - PDUItable$long.gp2
            PDUItable$short.gp2[PDUItable$short.gp2<0] <- 0
        }
        if(length(gp1)>1){
            PDUItable$total.mean.gp1 <- PDUItable$total.gp1
            PDUItable$total.mean.gp1[
                PDUItable$total.gp1<short_coverage_threshold] <- NA
            PDUItable$total.mean.gp1 <- rowMeans(PDUItable$total.mean.gp1, 
                                                 na.rm=TRUE)
            PDUItable$long.mean.gp1 <- PDUItable$long.gp1
            PDUItable$long.mean.gp1[
                PDUItable$long.gp1<long_coverage_threshold] <- NA
            PDUItable$long.mean.gp1 <- rowMeans(PDUItable$long.mean.gp1, 
                                                na.rm=TRUE)
            PDUItable$short.mean.gp1 <- rowMeans(PDUItable$short.gp1)
        }else{
            PDUItable$total.mean.gp1 <- PDUItable$total.gp1
            PDUItable$total.mean.gp1[
                PDUItable$total.gp1<short_coverage_threshold] <- NA
            PDUItable$long.mean.gp1 <- PDUItable$long.gp1
            PDUItable$long.mean.gp1[
                PDUItable$long.gp1<long_coverage_threshold] <- NA
            PDUItable$short.mean.gp1 <- PDUItable$short.gp1
        }
        if(is.character(gp2)){
            if(length(gp2)>1){
                PDUItable$total.mean.gp2 <- PDUItable$total.gp2
                PDUItable$total.mean.gp2[
                    PDUItable$total.gp2<short_coverage_threshold] <- NA
                PDUItable$total.mean.gp2 <- rowMeans(PDUItable$total.mean.gp2, 
                                                     na.rm=TRUE)
                PDUItable$long.mean.gp2 <- PDUItable$long.gp2
                PDUItable$long.mean.gp2[
                    PDUItable$long.gp2<long_coverage_threshold] <- NA
                PDUItable$long.mean.gp2 <- rowMeans(PDUItable$long.mean.gp2, 
                                                    na.rm=TRUE)
                PDUItable$short.mean.gp2 <- rowMeans(PDUItable$short.gp2)
            }else{
                PDUItable$total.mean.gp2 <- PDUItable$total.gp2
                PDUItable$total.mean.gp2[
                    PDUItable$total.gp2<short_coverage_threshold] <- NA
                PDUItable$long.mean.gp2 <- PDUItable$long.gp2
                PDUItable$long.mean.gp2[
                    PDUItable$long.gp2<long_coverage_threshold] <- NA
                PDUItable$short.mean.gp2 <- PDUItable$short.gp2
            }
        }
        
        data <- as.data.frame(mcols(PDUItable))
        if(is.character(gp2)){
            ##compensation
            ids <- data$short.mean.gp1==0 & data$short.mean.gp2==0
            ratio.gp1 <- data$long.mean.gp1/data$total.mean.gp1
            ratio.gp2 <- data$long.mean.gp2/data$total.mean.gp2
            ratio <- ifelse(ratio.gp1>ratio.gp2, ratio.gp1, ratio.gp2)
            short.mean.gp1 <- data$total.mean.gp1 * ratio - data$long.mean.gp1
            short.mean.gp2 <- data$total.mean.gp2 * ratio - data$long.mean.gp2
            short.mean.gp1[short.mean.gp1<0] <- 0
            short.mean.gp2[short.mean.gp2<0] <- 0
            PDUItable$short.mean.gp1[ids] <- short.mean.gp1[ids]
            PDUItable$short.mean.gp2[ids] <- short.mean.gp2[ids]
            PDUItable$total.mean.gp1[is.na(PDUItable$total.mean.gp1)] <- 0
            PDUItable$total.mean.gp2[is.na(PDUItable$total.mean.gp2)] <- 0
            PDUItable$long.mean.gp1[is.na(PDUItable$long.mean.gp1)] <- 0
            PDUItable$long.mean.gp2[is.na(PDUItable$long.mean.gp2)] <- 0
            PDUItable$short.mean.gp1[is.na(PDUItable$short.mean.gp1)] <- 0
            PDUItable$short.mean.gp2[is.na(PDUItable$short.mean.gp2)] <- 0
            PDUItable$test_status <- 
                PDUItable$total.mean.gp1>long_coverage_threshold & 
                PDUItable$total.mean.gp2>long_coverage_threshold & 
                PDUItable$long.mean.gp1>long_coverage_threshold & 
                PDUItable$long.mean.gp2>long_coverage_threshold
            data <- as.data.frame(mcols(PDUItable))
            pval <- apply(data[,c("long.mean.gp1", 
                                  "short.mean.gp1", 
                                  "long.mean.gp2", 
                                  "short.mean.gp2")], 
                          1, function(.d) 
                              fisher.test(matrix(floor(.d), ncol=2))$p.value)
            PDUI.gp1 <- 
                data$long.mean.gp1/(data$long.mean.gp1 + data$short.mean.gp1)
            PDUI.gp2 <- 
                data$long.mean.gp2/(data$long.mean.gp2 + data$short.mean.gp2)
            dPDUI <- PDUI.gp2 - PDUI.gp1
            PDUItable$PDUI.gp1 <- PDUI.gp1
            PDUItable$PDUI.gp2 <- PDUI.gp2
            PDUItable$dPDUI <- dPDUI
            PDUItable$logFC <- 
                abs(log2(PDUI.gp1+.Machine$double.xmin) - 
                        log2(PDUI.gp2+.Machine$double.xmin)) >= 
                PDUI_logFC_cutoff
        }else{
            PDUItable$test_status <- 
                PDUItable$total.mean.gp1>long_coverage_threshold & 
                PDUItable$long.mean.gp1>long_coverage_threshold
            data.long <- PDUItable$long.gp1
            data.short <- PDUItable$short.gp1
            data <- log2(cbind(data.long, data.short)+0.1)
            treatments <- cbind(long=c(rep(c(1,0),c(ncol(data.long), 
                                                    ncol(data.short)))), 
                                short=c(rep(c(0,1), c(ncol(data.long), 
                                                      ncol(data.short)))))
            design <- model.matrix(~-1+treatments)
            colnames(design) <- c("long", "short")
            fit <- lmFit(data, design)
            contrast.matrix<-makeContrasts(contrasts="long-short",
                                           levels=design)
            fit <- contrasts.fit(fit,contrast.matrix)
            fit <- eBayes(fit)
            data <- topTable(fit, number=nrow(fit), sort.by="none")
            pval <- data$P.Value
            PDUItable$dPDUI <- PDUItable$PDUI.gp1 <- 
                PDUItable$long.mean.gp1/(PDUItable$long.mean.gp1 + 
                                             PDUItable$short.mean.gp1)
            PDUItable$logFC <- TRUE
        }
    }
    PDUItable$pval <- pval
    PDUItable$adjPval <- p.adjust(pval, method="BH")
    
    PDUItable$filterPass <- PDUItable$adjPval<adjusted.P_val.cutoff & 
        abs(PDUItable$dPDUI)>dPDUI_cutoff & 
        PDUItable$test_status & 
        PDUItable$logFC
    PDUItable$source <- NULL
    PDUItable$logFC <- NULL
    
    PDUItable <- PDUItable[order(PDUItable$filterPass, 
                                 abs(PDUItable$dPDUI), 
                                 decreasing=TRUE)]
}