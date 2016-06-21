\name{distalAdj}
\alias{distalAdj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
adjust distal CP sites by cleanUpdTSeq
}
\description{
adjust distal CP sites by cleanUpdTSeq
}
\usage{
distalAdj(distalCPs, classifier, classifier_cutoff, shift_range, genome, step=1)
}
\arguments{
  \item{distalCPs}{
    the output of \link{searchDistalCPs}
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
  \item{genome}{
    a \link[BSgenome]{BSgenome} object
}
  \item{step}{
    adjust step, default 1, means adjust by each base by cleanUpdTSeq.
}
}
\value{
a list could be input of \code{\link{searchProximalCPs}}
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{searchDistalCPs}}, \code{\link{PAscore2}}
}
\keyword{ misc }