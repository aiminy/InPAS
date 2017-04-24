\name{proximalAdjByPWM}
\alias{proximalAdjByPWM}
\title{
    adjust the proximal CP sites by PWM
}
\description{
    adjust the proximal CP sites by polyA Position Weight Metrix. It only need 
    the PWM get match in upstream or downstream shift_range nr.
}
\usage{
proximalAdjByPWM(idx, PolyA_PWM, seqnames, starts, strands, genome, 
                 shift_range, search_point_START)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{idx}{
    the offset of positions of CP sites
}
  \item{PolyA_PWM}{
    polyA PWM
}
  \item{seqnames}{
    sequence names
}
  \item{starts}{
    start position in the genome
}
  \item{strands}{
    strands
}
  \item{genome}{
    an \link[BSgenome]{BSgenome} object
}
  \item{shift_range}{
    the shift range of PWM hits
}
  \item{search_point_START}{
    Not use
}
}
\value{
    the offset of positions of CP sites after filter
}
\details{
    the hits is searched by \link[Biostrings]{matchPWM} and the cutoff is 70\%
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{proximalAdjByCleanUpdTSeq}}, \code{\link{proximalAdj}},
    \code{\link{PAscore}}
}
\keyword{ misc }