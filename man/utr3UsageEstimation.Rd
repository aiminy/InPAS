\name{utr3UsageEstimation}
\alias{utr3UsageEstimation}
\title{
    estimation of 3'UTR usage for each region
}
\description{
    estimation of 3'UTR usage for short form and long form
}
\usage{
    utr3UsageEstimation(CPsites, coverage, genome, utr3, 
        gp1, gp2=NULL, 
        short_coverage_threshold = 10, 
        long_coverage_threshold = 2, 
        adjusted.P_val.cutoff = 0.05, 
        dPDUI_cutoff = 0.3, 
        PDUI_logFC_cutoff=0.59, BPPARAM=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{CPsites}{outputs of \link{CPsites}}
    \item{coverage}{coverage for each sample, outputs of 
        \link{coverageFromBedGraph}}
    \item{genome}{an object of \link[BSgenome]{BSgenome}}
    \item{utr3}{output of \link{utr3Annotation}}
    \item{gp1}{tag names involved in group 1}
    \item{gp2}{tag names involved in group 2}
    \item{short_coverage_threshold}{cutoff threshold for coverage in thre 
        region of short form}
    \item{long_coverage_threshold}{cutoff threshold for coverage in thre region 
        of long form}
    \item{adjusted.P_val.cutoff}{cutoff value for adjusted p.value}
    \item{dPDUI_cutoff}{cutoff value for differenctial 
        PAS(polyadenylation signal) usage index}
    \item{PDUI_logFC_cutoff}{cutoff value for log2 fold change of 
        PAS(polyadenylation signal) usage index}
    \item{BPPARAM}{An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.}
}
\value{
    return an object of GRanges
}
\author{
    Jianhong Ou
}
\examples{
    if(interactive()){
        library(BSgenome.Mmusculus.UCSC.mm10)
        path <- file.path(find.package("InPAS"), "extdata")
        bedgraphs <- file.path(path, "Baf3.extract.bedgraph")
        data(utr3.mm10)
        tags <- "Baf3"
        genome <- BSgenome.Mmusculus.UCSC.mm10
        coverage <- 
            coverageFromBedGraph(bedgraphs, tags, genome, hugeData=FALSE)
        CP <- CPsites(coverage=coverage, gp1=tags, gp2=NULL, genome=genome, 
            utr3=utr3.mm10, coverage_threshold=5, long_coverage_threshold=5)
        res <- utr3UsageEstimation(CP, coverage, 
            utr3.mm10, genome, gp1=tags, gp2=NULL)
    }
}
\keyword{ misc }