\name{UTR3usage}
\alias{UTR3usage}
\title{
    calculate the usage of long and short form of UTR3
}
\description{
    calculate the usage of long and short form of UTR3 for the results of \link{CPsites}
}
\usage{
UTR3usage(CPsites, coverage, hugeData, BPPARAM = NULL, phmm = FALSE)
}
\arguments{
  \item{CPsites}{
    outputs of \link{CPsites}
}
  \item{coverage}{
    coverage for each sample, outputs of 
        \link{coverageFromBedGraph}
}
  \item{hugeData}{
    is this dataset consume too much memory? 
        if it is TRUE, the coverage will be saved into tempfiles.
}
  \item{BPPARAM}{
    An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.
}
  \item{phmm}{
    prepare data for singleSample analysis? default is FALSE
}
}
\value{
    GRanges object
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{CPsites}}
}
\keyword{ misc }