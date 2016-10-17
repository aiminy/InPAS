\name{singleSampleAnalyze}
\alias{singleSampleAnalyze}
\title{
    do analysis for single sample
}
\description{
    do analysis for single sample by a hidden Markov model
}
\usage{
singleSampleAnalyze(UTR3eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UTR3eset}{
    must be the output of \code{\link{getUTR3eSet}}
}
}
\details{
    the test will be performed by a two states hidden Markov model.
}
\value{
a metrix of test results
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{UTR3eSet}}, \code{\link{getUTR3eSet}}, \code{\link{depmix}}
}
\examples{
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "eset.MAQC.rda"))
    res <- singleSampleAnalyze(eset)
}
\keyword{ misc }