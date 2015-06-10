\name{lastCDSusage}
\alias{lastCDSusage}
\title{
    extract coverage of last CDS exon region
}
\description{
    extract coverage of last CDS exon region
}
\usage{
lastCDSusage(CDS, coverage, hugeData, BPPARAM=NULL, phmm=FALSE)
}
\arguments{
  \item{CDS}{
    GRanges object of CDS
}
\item{coverage}{
    output of coverageFromBedGraph
}
  \item{hugeData}{
    is it a huge dataset?
}
  \item{BPPARAM}{
    An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.
}
  \item{phmm}{
    prepare data for singleSample analysis?
}
}
\value{
    the average coverage of last CDS for each transcript
}
\author{
    Jianhong Ou
}
\keyword{ misc }