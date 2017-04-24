\name{proximalAdj}
\alias{proximalAdj}
\title{
    adjust the proximal CP sites
}
\description{
    adjust the proximal CP sites by PolyA PWM and cleanUpdTSeq
}
\usage{
proximalAdj(CPs, MINSIZE, PolyA_PWM, genome, classifier, classifier_cutoff, 
            shift_range, search_point_START, step=1)
}

\arguments{
  \item{CPs}{
    the outputs of \link{searchProximalCPs}
}
  \item{MINSIZE}{
    min size for short from
}
  \item{PolyA_PWM}{
    PolyA position weight metrix
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
\value{
    keep same as \link{searchProximalCPs}, which can be handled by \link{polishCPs}.
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{searchProximalCPs}}, \code{\link{polishCPs}},
    \code{\link{proximalAdjByPWM}}, \code{\link{proximalAdjByCleanUpdTSeq}},
    \code{\link{PAscore}},\code{\link{PAscore2}}
}
\keyword{ misc }