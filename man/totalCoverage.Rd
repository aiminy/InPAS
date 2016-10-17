\name{totalCoverage}
\alias{totalCoverage}
\title{
    total coverage
}
\description{
    for huge dataset, it will read in the coverage from tmp files and merge them by groups
}
\usage{
totalCoverage(coverage, genome, hugeData, groupList=NULL)
}
\arguments{
  \item{coverage}{
    coverage for each sample, outputs of 
        \link{coverageFromBedGraph}
}
  \item{genome}{
    an object of \link[BSgenome]{BSgenome}
}
  \item{hugeData}{
    hugeData or not
}
  \item{groupList}{
    tag names involved in each groups
}
}
\value{
    a coverage list
}
\author{
    Jianhong Ou
}
\keyword{ misc }