\name{PAscore}
\alias{PAscore}
\title{
calculate the CP score
}
\description{
calculate the CP score by PWM
}
\usage{
PAscore(seqname, pos, str, idx, PWM, genome, ups = 50, dws = 50)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{PWM}{
    polyA position weight matrix
}
  \item{genome}{
    an object of \link[BSgenome]{BSgenome}
}
  \item{ups}{
    upstream base
}
  \item{dws}{
    downstream base
}
}
\value{
    idx list after filter
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{PAscore2}}
}
\keyword{ misc }