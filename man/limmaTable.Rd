\name{limmaTable}
\alias{limmaTable}
\title{
    use limma to analyze the PDUI
}
\description{
    use limma to analyze the PDUI
}
\usage{
limmaTable(fit, UTR3eset, coef=1, ...)
}
\arguments{
  \item{fit}{
    fit should be an object of class \link[limma]{MArrayLM-class} as produced by \link[limma]{lmFit} and \link[limma]{eBayes}.
}
  \item{UTR3eset}{
    an \link{UTR3eSet} object
}
  \item{coef}{
    column number or column name specifying which coefficient or contrast of the linear model is of interest. see details at \link[limma]{topTable}
}
  \item{...}{
    any parameters can be passed to \link[limma]{topTable}
}
}
\value{
    a dataframe output by \link[limma]{topTable}. The output is reorder by the transcript name of inputs.
}
\author{
    Jianhong Ou
}

\seealso{
    \code{\link{singleSampleAnalyze}},\code{\link{singleGroupAnalyze}},
    \code{\link{fisher.exact.test}}
}
\examples{
    library(limma)
    path <- file.path(find.package("InPAS"), "extdata")
    load(file.path(path, "eset.MAQC.rda"))
    tags <- colnames(eset$PDUI.log2)
    g <- factor(gsub("\\\\..*$", "", tags))
    design <- model.matrix(~-1+g)
    colnames(design) <- c("Brain", "UHR")
    contrast.matrix <- makeContrasts(contrasts="Brain-UHR",levels=design)
    fit <- limmaAnalyze(eset, design, contrast.matrix)
    limmaTable(fit, eset)
}
\keyword{ misc }
