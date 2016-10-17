\name{optimalSegmentation}
\alias{optimalSegmentation}
\title{
    calculate SSE
}
\description{
    calculate SSE values
}
\usage{
optimalSegmentation(.ele, search_point_START, search_point_END, n = 1, savedID = NA)
}
\arguments{
  \item{.ele}{
    3UTR coverage
}
  \item{search_point_START}{
    start position to calculate
}
  \item{search_point_END}{
    end position to calculate
}
  \item{n}{
    the length of output
}
  \item{savedID}{
    the proximal CPsites for noval distal events
}
}
\value{
a list of SSE and idx
}
\author{
    Jianhong Ou
}
\keyword{ misc }