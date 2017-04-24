\name{UTR3eSet-class}
\docType{class}
\alias{UTR3eSet-class}
\alias{UTR3eSet}
\alias{$,UTR3eSet-method}
\alias{$<-,UTR3eSet-method}

\title{
    Class \code{UTR3eSet}
}
\description{
    An object of class \code{UTR3eSet} represents the results of 3UTR usage
}
\section{Objects from the Class}{
    Objects can be created by calls of the form
    \code{new("UTR3eSet", usage, PDUI, PDUI.log2, short, long, signals, testRes)}
}
\section{Slots}{
  \describe{
    \item{\code{usage}}{an \link[GenomicRanges]{GRanges} object with CP sites info.}
    
    \item{\code{PDUI}}{a matrix of PDUI}
    
    \item{\code{PDUI.log2}}{log2 transformed PDUI matrix}
    
    \item{\code{short}}{a matrix of usage of short form}
    
    \item{\code{long}}{a matrix of usage of long form}
    
    \item{\code{signals}}{signals used for single sample}
    
    \item{\code{testRes}}{a matrix of test results of \link{testUsage}}
  }
}
\section{Methods}{
  \describe{
    \item{$, $<-}{Get or set the slot of \link{UTR3eSet}}
    \item{as("UTR3eSet", "ExpressionSet")}{Convert a UTR3eSet to an \link[Biobase]{ExpressionSet}.}
    \item{as("UTR3eSet", "GRanges")}{Convert a UTR3eSet to an \link[GenomicRanges]{GRanges}.}
  }
}
\author{
    Jianhong Ou
}
\keyword{ classes }