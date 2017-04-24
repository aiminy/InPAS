coverageCDS <- function(coverage, CDS, genome, gp1, gp2, hugeData){
    totalCov <- totalCoverage(coverage, genome, hugeData, groupList=list(gp1=gp1, gp2=gp2))
    cdsTotalCov <- UTR3TotalCoverage(CDS, totalCov)
    cdsTotalCov <- lapply(cdsTotalCov, function(.ele){
        lapply(.ele, function(.e){
            .e <- colSums(.e)/nrow(.e)
        })
    })
    cdsTotalCov <- do.call(rbind, lapply(cdsTotalCov, do.call, what=rbind))
    if(!hugeData){
        gp1 <- cdsTotalCov[, gp1, drop=FALSE]
        gp2 <- cdsTotalCov[, gp2, drop=FALSE]
        gp1 <- rowMeans(gp1)
        gp2 <- rowMeans(gp2)
        cdsTotalCov <- cbind(gp1=gp1, gp2=gp2)
    }
    cdsTotalCov
}

utr3UsageEstimation <- function(CPsites, coverage, genome, utr3, 
                                gp1, gp2=NULL, 
                                short_coverage_threshold=10, 
                                long_coverage_threshold=2, 
                                adjusted.P_val.cutoff=0.05, 
                                dPDUI_cutoff=0.3, 
                                PDUI_logFC_cutoff=0.59, 
                                BPPARAM=NULL){
    if(!all(c(gp1,gp2) %in% names(coverage))) 
        stop("gp1 and gp2 must be in names of coverage")
    if(missing(coverage) || missing(CPsites))
        stop("CPsites and coverage is required.")
    if(length(gp2)>0){
        if(missing(utr3) || missing(genome)){
            stop("utr3 and genome is required.")
        }
        if(class(genome)!="BSgenome")
            stop("genome must be an object of BSgenome.")
        if(class(utr3)!="GRanges" | 
               !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))){
            stop("utr3 must be output of function of utr3Annotation")
        }
    }
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
            ##compensation by last CDS
            ids <- data$short.mean.gp1<=0 & data$short.mean.gp2<=0
            if(sum(ids)>0){
                ##get coverage of last CDS
                CDS <- utr3[utr3$feature=="CDS"]
                CDS <- CDS[CDS$transcript %in% unique(CPsites$transcript)]
                cds.id <- match(data$transcript[ids], CDS$transcript)
                cds.id.na <- is.na(cds.id)
                CDS <- CDS[cds.id[!cds.id.na]]
                if(length(CDS)>0){
                    ##use last CDS coverage replace the short.mean.gp1 and gp2
                    CDS.cov <- coverageCDS(coverage, CDS, genome, 
                                           gp1, gp2, hugeData)
                    CDS.cov <- CDS.cov[match(names(CDS), rownames(CDS.cov)), , 
                                       drop=FALSE]
                    rownames(CDS.cov) <- names(CDS)
                    CDS.cov[is.na(CDS.cov)] <- 0
                    PDUItable$short.mean.gp1[ids][!cds.id.na] <- CDS.cov[,"gp1"] - 
                        PDUItable$long.mean.gp1[ids][!cds.id.na]
                    PDUItable$short.mean.gp2[ids][!cds.id.na] <- CDS.cov[,"gp2"] - 
                        PDUItable$long.mean.gp2[ids][!cds.id.na]
                    PDUItable$short.mean.gp1[
                        is.na(PDUItable$short.mean.gp1) | 
                            PDUItable$short.mean.gp1<0] <- 0
                    PDUItable$short.mean.gp2[
                        is.na(PDUItable$short.mean.gp2) | 
                            PDUItable$short.mean.gp2<0] <- 0
                }
            }
            
            ##compensation by ratio
            #data <- as.data.frame(mcols(PDUItable))
            #ids <- data$short.mean.gp1<=0 & data$short.mean.gp2<=0
            #ratio.gp1 <- data$long.mean.gp1/data$total.mean.gp1
            #ratio.gp2 <- data$long.mean.gp2/data$total.mean.gp2
            #ratio <- ifelse(ratio.gp1>ratio.gp2, ratio.gp1, ratio.gp2)
            #short.mean.gp1 <- data$total.mean.gp1 * ratio - data$long.mean.gp1
            #short.mean.gp2 <- data$total.mean.gp2 * ratio - data$long.mean.gp2
            #short.mean.gp1[short.mean.gp1<0] <- 0
            #short.mean.gp2[short.mean.gp2<0] <- 0
            #PDUItable$short.mean.gp1[ids] <- short.mean.gp1[ids]
            #PDUItable$short.mean.gp2[ids] <- short.mean.gp2[ids]
            PDUItable$total.mean.gp1[is.na(PDUItable$total.mean.gp1)] <- 0
            PDUItable$total.mean.gp2[is.na(PDUItable$total.mean.gp2)] <- 0
            PDUItable$long.mean.gp1[is.na(PDUItable$long.mean.gp1)] <- 0
            PDUItable$long.mean.gp2[is.na(PDUItable$long.mean.gp2)] <- 0
            PDUItable$short.mean.gp1[is.na(PDUItable$short.mean.gp1)] <- 0
            PDUItable$short.mean.gp2[is.na(PDUItable$short.mean.gp2)] <- 0
            PDUItable$test_status <- 
                (PDUItable$total.mean.gp1>=long_coverage_threshold | 
                     PDUItable$total.mean.gp2>=long_coverage_threshold) & 
                (PDUItable$long.mean.gp1>=long_coverage_threshold | 
                     PDUItable$long.mean.gp2>=long_coverage_threshold) &
                !(PDUItable$total.mean.gp1<long_coverage_threshold & 
                      PDUItable$long.mean.gp1<long_coverage_threshold)
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