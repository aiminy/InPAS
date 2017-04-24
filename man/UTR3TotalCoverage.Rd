\name{UTR3TotalCoverage}
\alias{UTR3TotalCoverage}
\title{
    extract coverage of 3UTR for CP sites prediction
}
\description{
    extract 3UTR coverage from totalCov according and GRanges object utr3.
}
\usage{
UTR3TotalCoverage(utr3, totalCov, gcCompensation = NA, 
                  mappabilityCompensation = NA, 
                  FFT = FALSE, fft.sm.power = 20)
}
\arguments{
  \item{utr3}{
    an GRanges object. must be the output of \link{utr3Annotation}
}
  \item{totalCov}{
    total coverage of each sample. must be the output of \link{totalCoverage}
}
  \item{gcCompensation}{
  GC compensation vector. Not support yet.
}
  \item{mappabilityCompensation}{
    mappability compensation vector. Not support yet.
}
  \item{FFT}{
    Use FFT smooth or not.
}
  \item{fft.sm.power}{
    the cut-off frequence of FFT smooth.
}
}

\value{
    a list. level 1: chromosome; level 2: each transcripts; level3: data matrix
}
\author{
    Jianhong Ou
}
\keyword{ misc }