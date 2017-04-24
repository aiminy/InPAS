\name{zScoreThrethold}
\alias{zScoreThrethold}

\title{
calculate local background cutoff value
}
\description{
calculate local background cutoff value based on z-score
}
\usage{
zScoreThrethold(background, introns, totalCov, utr3, z = 2)
}

\arguments{
  \item{background}{
    background range
}
  \item{introns}{
    GRanges of introns
}
  \item{totalCov}{
    total coverage of output of \link{totalCoverage}
}
  \item{utr3}{
    output of \link{utr3Annotation}
}
  \item{z}{
    z score cut off value
}
}
\value{
a numeric vector
}
\author{
Jianhong Ou
}
\keyword{ misc }
