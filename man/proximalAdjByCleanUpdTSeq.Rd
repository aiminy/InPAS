\name{proximalAdjByCleanUpdTSeq}
\alias{proximalAdjByCleanUpdTSeq}
\title{
    adjust the proximal CP sites by cleanUpdTseq
}
\description{
    adjust the proximal CP sites by cleanUpdTseq
}
\usage{
proximalAdjByCleanUpdTSeq(idx.list, cov_diff.list, seqnames, starts, strands, 
                        genome, classifier, classifier_cutoff, 
                        shift_range, search_point_START, step=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{idx.list}{
    the offset of positions of CP sites
}
  \item{cov_diff.list}{
    the SSE values
}
  \item{seqnames}{
    sequence names
}
  \item{starts}{
    starts
}
  \item{strands}{
    strands
}
  \item{genome}{
    a \link[BSgenome]{BSgenome} object
}
  \item{classifier}{
    cleanUpdTSeq classifier
}
  \item{classifier_cutoff}{
    cutoff value of the classiffier
}
  \item{shift_range}{
    the searching range for the better CP sites
}
  \item{search_point_START}{
    just in case there is no better CP sites
}
  \item{step}{
    adjust step, default 1, means adjust by each base by cleanUpdTSeq.
}
}
\details{
    the step for calculating is 10, can not do every base base it is really very slow.
}
\value{
    the offset of positions of CP sites after filter
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{proximalAdjByPWM}}, \code{\link{proximalAdj}},
    \code{\link{PAscore2}}
}
\keyword{ misc }