\name{searchDistalCPs}
\alias{searchDistalCPs}
\title{
    search distal CP sites
}
\description{
    search distal CP sites
}
\usage{
searchDistalCPs(chr.cov.merge, conn_next_utr3, 
            curr_UTR, window_size, 
            depth.weight, 
            long_coverage_threshold, 
            background, z2s)
}
\arguments{
  \item{chr.cov.merge}{
    coverage of current chromosome
}
  \item{conn_next_utr3}{
    joint to next 3UTR or not (used for \link{removeUTR3__UTR3})
}
  \item{curr_UTR}{
    GRanges of current 3UTR
}
  \item{window_size}{
    window size
}
  \item{depth.weight}{
    output of \link{depthWeight}
}
  \item{long_coverage_threshold}{
    cutoff value for coverage of long form 3UTR
}
  \item{background}{
    local background range
}
  \item{z2s}{
    cut off background scores. see \link{zScoreThrethold}
}
}
\value{
a list
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{distalAdj}}, \code{\link{PAscore2}}
}
\keyword{ misc }