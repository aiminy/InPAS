\name{singleGroupAnalyze}
\alias{singleGroupAnalyze}
\title{
    do analysis for single group samples
}
\description{
    do analysis for single group samples by anova test
}
\usage{
singleGroupAnalyze(UTR3eset)
}
\arguments{
  \item{UTR3eset}{
    must be the output of \code{\link{getUTR3eSet}}
}
}
\value{
    a metrix of test results
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{UTR3eSet}}, \code{\link{getUTR3eSet}}
}
\examples{
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "eset.MAQC.rda"))
    res <- singleGroupAnalyze(eset)
}
\keyword{ misc }