\name{testUsage}
\alias{testUsage}
\title{
 do test for dPDUI
}
\description{
    do test for dPDUI
}
\usage{
testUsage(CPsites, coverage, genome, utr3, BPPARAM=NULL, 
          method=c("limma", "fisher.exact", 
                    "singleSample", "singleGroup"),
          normalize=c("none", "quantiles", "quantiles.robust",
                        "mean", "median"),
          design, contrast.matrix, coef=1,
          gp1, gp2)
}
\arguments{
    \item{CPsites}{outputs of \link{CPsites}}
    \item{coverage}{coverage for each sample, outputs of 
        \link{coverageFromBedGraph}}
    \item{genome}{an object of \link[BSgenome]{BSgenome}}
    \item{utr3}{output of \link{utr3Annotation}}
  \item{BPPARAM}{
    An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.
}
  \item{method}{
    test method. see \code{\link{singleSampleAnalyze}},\code{\link{singleGroupAnalyze}},
    \code{\link{fisher.exact.test}},\code{\link{limmaAnalyze}} 
}
  \item{normalize}{
    normalization method
}
  \item{design}{
    the design matrix of the experiment, with rows corresponding to arrays and columns to coefficients to be estimated. Defaults to the unit vector meaning that the arrays are treated as replicates. see \link[stats]{model.matrix}
}
  \item{contrast.matrix}{
    numeric matrix with rows corresponding to coefficients in fit and columns containing contrasts. May be a vector if there is only one contrast. see \link[limma]{makeContrasts}
}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest. see more \link[limma]{topTable}. default value: 1}
  \item{gp1}{
    tag names involved in group 1
}
  \item{gp2}{
    tag names involved in group 2
}
}
\details{
    if method is "limma", design matrix and contrast is required.
    if method is "fisher.exact", gp1 and gp2 is required.
}
\value{
    a list with test results. the output of test results is a matrix.
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{singleSampleAnalyze}},\code{\link{singleGroupAnalyze}},
    \code{\link{fisher.exact.test}},\code{\link{limmaAnalyze}} 
}
\examples{
    library(limma)
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "CPs.MAQC.rda"))
    load(file.path(path, "coverage.MAQC.rda"))
    library(BSgenome.Hsapiens.UCSC.hg19)
    data(utr3.hg19)
    tags <- names(coverage)
    g <- factor(gsub("\\\\..*$", "", tags))
    design <- model.matrix(~-1+g)
    colnames(design) <- c("Brain", "UHR")
    contrast.matrix<-makeContrasts(contrasts="Brain-UHR",levels=design)
    res <- testUsage(CPsites=CPs, 
                 coverage=coverage, 
                 genome=BSgenome.Hsapiens.UCSC.hg19,
                 utr3=utr3.hg19, 
                 method="limma",
                 design=design,
                 contrast.matrix=contrast.matrix)
}
\keyword{ misc }