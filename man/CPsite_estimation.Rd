\name{CPsite_estimation}
\alias{CPsite_estimation}
\title{
    estimate the cpsites
}
\description{
    estimate the cpsites for a giving chromosome
}
\usage{
CPsite_estimation(chr.cov, utr3, MINSIZE, window_size, search_point_START, 
search_point_END, cutStart, cutEnd, adjust_distal_polyA_end, 
background, z2s, coverage_threshold, long_coverage_threshold, 
PolyA_PWM, classifier, classifier_cutoff, shift_range, 
depth.weight, genome, step=1,
tmpfolder=NULL, silence=TRUE)
}
\arguments{
  \item{chr.cov}{
    coverage list for one chromosome
}
  \item{utr3}{
    output of utr3Annotaion
}
  \item{MINSIZE}{
    min size of short form
}
  \item{window_size}{
    window size
}
  \item{search_point_START}{
    search start point
}
  \item{search_point_END}{
    search end point
}
  \item{cutStart}{
    cut from start
}
  \item{cutEnd}{
    cut from end
}
  \item{adjust_distal_polyA_end}{
    adjust distal site or not
}
  \item{background}{
    how to get the local background
}
  \item{z2s}{
    output of \link{zScoreThrethold}
}
  \item{coverage_threshold}{
    cutoff value for coverage
}
  \item{long_coverage_threshold}{
    cutoff value for long form
}
  \item{PolyA_PWM}{
    polyA PWM
}
  \item{classifier}{
    classifier
}
  \item{classifier_cutoff}{
    classifier cutoff
}
  \item{shift_range}{
    shift range
}
  \item{depth.weight}{
    output of \link{depthWeight}
}
  \item{genome}{
    a \link[BSgenome]{BSgenome} object
}
  \item{step}{
    adjust step, default 1, means adjust by each base by cleanUpdTSeq.
}
  \item{tmpfolder}{
    temp folder could save and reload the analysis data for resume analysis.
  }
  \item{silence}{report progress or not. default not report.}
}
\value{
a data.frame
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{CPsites}}, \code{\link{searchProximalCPs}}, 
    \code{\link{proximalAdj}},
    \code{\link{proximalAdjByPWM}}, \code{\link{proximalAdjByCleanUpdTSeq}},
    \code{\link{PAscore}},\code{\link{PAscore2}}
}
\keyword{ misc }