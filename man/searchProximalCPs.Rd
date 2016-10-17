\name{searchProximalCPs}
\alias{searchProximalCPs}
\title{
search proximal CPsites
}
\description{
    search proximal CPsites
}
\usage{
searchProximalCPs(CPs, curr_UTR, window_size, 
                MINSIZE, cutEnd, 
                search_point_START, 
                search_point_END,
                two_way=FALSE)
}
\arguments{
  \item{CPs}{
    output of \link{searchDistalCPs} or \link{distalAdj}
}
  \item{curr_UTR}{
    GRanges of current 3UTR
}
  \item{window_size}{
    window size
}
  \item{MINSIZE}{
    MINSIZE for short form
}
  \item{cutEnd}{
    how many nucleotides should be removed  from the end before 
        search, 0.1 means 10 percent.
}
  \item{search_point_START}{
    start point for searching
}
  \item{search_point_END}{
    end point for searching
}
  \item{two_way}{
    Search the proximal site from both direction or not.
}
}

\value{
a list
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{proximalAdj}}, \code{\link{polishCPs}},
    \code{\link{proximalAdjByPWM}}, \code{\link{proximalAdjByCleanUpdTSeq}},
    \code{\link{PAscore}},\code{\link{PAscore2}}
}
\keyword{ misc }