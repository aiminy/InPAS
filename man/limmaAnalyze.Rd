\name{limmaAnalyze}
\alias{limmaAnalyze}
\title{
    use limma to analyze the PDUI
}
\description{
    use limma to analyze the PDUI
}
\usage{
limmaAnalyze(UTR3eset, design, contrast.matrix, coef=1, robust=FALSE, ...)
}
\arguments{
  \item{UTR3eset}{
    an \link{UTR3eSet} object
}
  \item{design}{
    the design matrix of the experiment, with rows corresponding to arrays and columns to coefficients to be estimated. Defaults to the unit vector meaning that the arrays are treated as replicates. see \link[stats]{model.matrix}
}
  \item{contrast.matrix}{
    numeric matrix with rows corresponding to coefficients in fit and columns containing contrasts. May be a vector if there is only one contrast. see \link[limma]{makeContrasts}
}
    \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest. see more \link[limma]{topTable}. default value: 1}
    \item{robust}{logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances?}
    \item{...}{other arguments are passed to lmFit.}
}
\value{
    fit results of eBayes by limma. It is an object of class MArrayLM containing everything found in fit. see \link[limma]{eBayes}
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
    res <- limmaAnalyze(eset, design, contrast.matrix)
    head(res)
}
\keyword{ misc }
