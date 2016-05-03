\name{getCov}
\alias{getCov}
\title{
    extract coverage from bedgraph file
}
\description{
    extract coverage from bedgraph file
}
\usage{
getCov(bedgraph, genome, seqLen)
}
\arguments{
  \item{bedgraph}{
    bedGraph file names
}
  \item{genome}{
    an object \link[BSgenome]{BSgenome}
}
  \item{seqLen}{
    lengthes of each chromosome
}
}
\value{
    a Rle object for a sample coverage 
}
\author{
    Jianhong Ou
}

\seealso{
    \code{\link{coverageFromBedGraph}}
}
\keyword{ misc }
