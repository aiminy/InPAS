\name{usage4plot}
\alias{usage4plot}
\title{
    prepare coverage data and fitting data for plot
}
\description{
    prepare coverage data and fitting data for plot
}
\usage{
    usage4plot(gr, coverage, proximalSites, genome, groupList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gr}{an object of GRanges}
    \item{coverage}{coverage for each sample}
    \item{proximalSites}{proximal sites}
    \item{genome}{an object of \link[BSgenome]{BSgenome}}
    \item{groupList}{the list of sample names}
}
\value{
    Formal class 'GRanges' [package "GenomicRanges"] with metadata: 
    \item{dat}{matrix, first column is the fit data, the other columns are 
        coverage data for each sample}
    \item{offset}{offset from the start of 3UTR}
}
\author{
    Jianhong Ou
}
\examples{
    library(BSgenome.Mmusculus.UCSC.mm10)
    path <- file.path(find.package("InPAS"), "extdata")
    bedgraphs <- c(file.path(path, "Baf3.extract.bedgraph"), 
               file.path(path, "UM15.extract.bedgraph"))
    coverage <- coverageFromBedGraph(bedgraphs, tags=c("Baf3", "UM15"), 
                                genome=Mmusculus, hugeData=FALSE)
    gr <- GRanges("chr6", IRanges(128846245, 128850081), strand="-")
    dat <- usage4plot(gr, coverage, proximalSites=128849148, Mmusculus)
    data <- dat$dat[[1]]
    op <- par(mfrow=c(3, 1))
    plot(data[,1], type="l", xlab="", ylab="The fitted value")
    abline(v=dat$offset)
    plot(data[,2], type="l", xlab="", ylab="Baf3")
    plot(data[,3], type="l", xlab="", ylab="UM15")
    par(op)
}
\keyword{ misc }
