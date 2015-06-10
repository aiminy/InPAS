\name{get.regions.coverage}
\alias{get.regions.coverage}
\title{
    claculate coverage for giving region
}
\description{
    claculate coverage for giving region
}
\usage{
get.regions.coverage(chr, utr3.regions.chr, 
                                 hugeData, coverage, phmm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
    chromosome
}
  \item{utr3.regions.chr}{
    the GRanges of region to be extracted
}
  \item{hugeData}{
    is it a huge dataset?
}
  \item{coverage}{
    output of coverageFromBedGraph
}
  \item{phmm}{
    prepare data for singleSample analysis?
}
}
\value{
GRanges with coverage data
}
\author{
    Jianhong Ou
}

\keyword{ misc }