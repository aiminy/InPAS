\name{depthWeight}
\alias{depthWeight}
\title{
    calculate the depth weight for each example
}
\description{
    calculate the depth weight for each example
}
\usage{
depthWeight(coverage, hugeData, groupList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coverage}{
    a list. output of \link{coverageFromBedGraph}
}
  \item{hugeData}{
    is it a huge dataset?
}
  \item{groupList}{
    group list for huge dataset
}
}
\value{
    a numeric vector with depth weight
}
\author{
    Jianhong Ou
}
\keyword{ misc }