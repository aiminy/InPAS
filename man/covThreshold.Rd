\name{covThreshold}
\alias{covThreshold}
\title{
    calculate the cutoff threshold of coverage 
}
\description{
    calculate the cutoff threshold of coverage for long form and short form
}
\usage{
covThreshold(coverage, genome, txdb, utr3, 
                         chr="chr1", hugeData, groupList)
}
\arguments{
  \item{coverage}{
    coverage for each sample, 
        output of \link{coverageFromBedGraph}
}
  \item{genome}{
    an object of \link[BSgenome]{BSgenome}
}
  \item{txdb}{
    an object of \link[GenomicFeatures]{TxDb}
}
  \item{utr3}{
    output of \link{utr3Annotation}
}
  \item{chr}{
    chromosome to be used for calculation, default is "chr1"
}
  \item{hugeData}{
    is this dataset consume too much memory? 
        if it is TRUE, the coverage will be saved into tempfiles.
}
  \item{groupList}{
    group list of tag names 
}
}
\value{
 a numeric vector
}
\author{
    Jianhong Ou
}

\seealso{
    \code{\link{CPsite_estimation}}
}
\keyword{ misc }