\name{prepare4GSEA}
\alias{prepare4GSEA}
\title{
prepare the files for GSEA analysis
}
\description{
output the log2 transformed delta PDUI txt file and chip file for GSEA analysis
}
\usage{
prepare4GSEA(eset, groupList, Preranked=TRUE,
             folder=".", 
             rnkFilename="InPAS.rnk",
             chipFilename="InPAS.chip", 
             dataFilename="dPDUI.txt",
             PhenFilename="group.cls")
}
\arguments{
  \item{eset}{
    a \link{UTR3eSet} object
}
  \item{groupList}{
  group list of tag names
}
  \item{Preranked}{
  logical value, out preranked or not
}
  \item{folder}{
    output folder
}
  \item{rnkFilename}{
    filename of preranked file
}
  \item{chipFilename}{
    filename of chip
}
  \item{dataFilename}{
    filename of dataset
}
  \item{PhenFilename}{
  filename of Phenotype labels
}
}
\value{
None
}
\author{
Jianhong Ou
}
\examples{
    file <- system.file("extdata", "eset.MAQC.rda", package="InPAS")
    load(file)
    gp1=c("Brain.auto", "Brain.phiX")
    gp2=c("UHR.auto", "UHR.phiX")
    groupList <- list(Brain=gp1, UHR=gp2)
    prepare4GSEA(eset, groupList=groupList, Preranked=FALSE)
}
\keyword{ misc }
