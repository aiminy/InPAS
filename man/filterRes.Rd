\name{filterRes}
\alias{filterRes}
\title{
 filter results
}
\description{
 filter results of \link{testUsage}
}
\usage{
filterRes(res,
          gp1, gp2, 
          background_coverage_threshold=2, 
          P.Value_cutoff=0.05,
          adj.P.Val_cutoff=0.05, 
          dPDUI_cutoff=0.3, 
          PDUI_logFC_cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
    output of \link{testUsage}
}
  \item{gp1}{tag names involved in group 1}
  \item{gp2}{tag names involved in group 2}
  \item{background_coverage_threshold}{
    background coverage cut off value.
    for each group, more than half of the long form should greater than background_coverage_threshold.
    for both group, at leat in one group, more than half of the short form should greater than background_coverage_threshold.
  }
  \item{P.Value_cutoff}{
    cutoff of P value
  }
  \item{adj.P.Val_cutoff}{
    cutoff of adjust P value
  }
  \item{dPDUI_cutoff}{
    cutoff of dPDUI
  }
  \item{PDUI_logFC_cutoff}{
    cutoff of PDUI log2 transformed fold change
  }
}
\value{
a data.frame
}
\author{
    Jianhong Ou
}
\seealso{
    \code{\link{testUsage}}
}
\examples{
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "CPs.MAQC.rda"))
    load(file.path(path, "coverage.MAQC.rda"))
    library(BSgenome.Hsapiens.UCSC.hg19)
    data(utr3.hg19)
    res <- testUsage(CPsites=CPs, 
                  coverage=coverage, 
                  genome=BSgenome.Hsapiens.UCSC.hg19,
                  utr3=utr3.hg19, 
                  method="fisher.exact",
                  gp1=c("Brain.auto", "Brain.phiX"),
                  gp2=c("UHR.auto", "UHR.phiX"))
    filterRes(res, 
              gp1=c("Brain.auto", "Brain.phiX"),
              gp2=c("UHR.auto", "UHR.phiX"),
              background_coverage_threshold=2, 
              P.Value_cutoff=0.05,
              adj.P.Val_cutoff=0.05, 
              dPDUI_cutoff=0.3, 
              PDUI_logFC_cutoff=.59)
}
\keyword{ misc }
