\name{polishCPs}
\alias{polishCPs}
\title{
    polish the searching results of CP sites
}
\description{
    remove the multiple positions of CP sites for same 3UTRs and only keep the best
    CP sites for proximal and distal.
}
\usage{
polishCPs(CPs)
}
\arguments{
  \item{CPs}{
    output of \link{searchProximalCPs} or \link{proximalAdj}
}
}
\value{
    a matrix with columes: "fit_value", "Predicted_Proximal_APA", 
    "Predicted_Distal_APA", "utr3start", "utr3end", "type"
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{CPsite_estimation}}, \code{\link{searchProximalCPs}}, 
    \code{\link{proximalAdj}},
    \code{\link{proximalAdjByPWM}}, \code{\link{proximalAdjByCleanUpdTSeq}},
    \code{\link{PAscore}},\code{\link{PAscore2}}
}
\keyword{ misc }