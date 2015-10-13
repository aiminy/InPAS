\name{fisher.exact.test}
\alias{fisher.exact.test}
\title{
    do fisher exact test for two group datasets
}
\description{
    do fisher exact test for two group datasets
}
\usage{
fisher.exact.test(UTR3eset, gp1, gp2)
}
\arguments{
  \item{UTR3eset}{
    output of \link{getUTR3eSet}
}
  \item{gp1}{
    tag names of group 1
}
  \item{gp2}{
    tag names of group 2
}
}
\value{
a matrix of test results
}
\author{
    Jianhong Ou
}

\seealso{
    \code{\link{singleSampleAnalyze}},\code{\link{singleGroupAnalyze}},
    \code{\link{limmaAnalyze}} 
}
\examples{
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "eset.MAQC.rda"))
    tags <- colnames(eset$PDUI.log2)    
    res <- fisher.exact.test(eset, gp1=tags[1:2], gp2=tags[3:4])
}
\keyword{ misc }