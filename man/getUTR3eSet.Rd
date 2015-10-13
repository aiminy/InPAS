\name{getUTR3eSet}
\alias{getUTR3eSet}
\title{
    prepare dataset for test
}
\description{
    Generate a UTR3eSet object with PDUI infomation for statistic test
}
\usage{
getUTR3eSet(CPsites, coverage, genome, utr3, 
            normalize=c("none", "quantiles", "quantiles.robust",
                      "mean", "median"),
            ...,
            BPPARAM=NULL, singleSample=FALSE)
}
\arguments{
  \item{CPsites}{
    outputs of \link{CPsites}
}
  \item{coverage}{
    coverage for each sample, outputs of 
        \link{coverageFromBedGraph}
}
  \item{genome}{
    an object of \link[BSgenome]{BSgenome}
}
  \item{utr3}{
    output of \link{utr3Annotation}
}
  \item{normalize}{
    normalization method
}
  \item{...}{
    parameter can be passed into \link[preprocessCore]{normalize.quantiles.robust}
}
  \item{BPPARAM}{
    An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.
}
  \item{singleSample}{
    prepare data for singleSample analysis? default is FALSE
}
}
\value{
    An object of \link{UTR3eSet} which contains following elements:
    
    usage: an GRanges object with CP sites info.
    
    PDUI: a matrix of PDUI
    
    PDUI.log2: log2 transformed PDUI matrix
    
    short: a matrix of usage of short form
    
    long: a matrix of usage of long form
    
    if singleSample is TRUE, one more element, signals, will be included.
}
\author{
    Jianhong Ou
}
\examples{
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "CPs.MAQC.rda"))
    load(file.path(path, "coverage.MAQC.rda"))
    library(BSgenome.Hsapiens.UCSC.hg19)
    data(utr3.hg19)
    getUTR3eSet(CPsites=CPs, 
                 coverage=coverage, 
                 genome=BSgenome.Hsapiens.UCSC.hg19,
                 utr3=utr3.hg19)
}

\keyword{ misc }