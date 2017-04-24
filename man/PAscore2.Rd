\name{PAscore2}
\alias{PAscore2}
\title{
 calculate the CP score
}
\description{
    calculate CP score by cleanUpdTSeq
}
\usage{
PAscore2(seqname, pos, str, idx, idx.gp, genome, classifier, classifier_cutoff)
}

\arguments{
  \item{seqname}{
    sequence names
}
  \item{pos}{
    genomic positions
}
  \item{str}{
    strands
}
  \item{idx}{
    offset postion
}
  \item{idx.gp}{
    group number of the offset postion
}
  \item{genome}{
    an object of \link[BSgenome]{BSgenome}
}
  \item{classifier}{
    a cleanUpdTSeq classifier
}
  \item{classifier_cutoff}{
    classifier cutoff value
}
}
\value{
a data.frame 
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{PAscore}}
}
\keyword{ misc }