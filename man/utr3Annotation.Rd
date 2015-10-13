\name{utr3Annotation}
\alias{utr3Annotation}
\title{
    extract 3'UTR from \link[GenomicFeatures]{TxDb} object
}
\description{
    extract 3'UTR from a \link[GenomicFeatures]{TxDb} object. The 3'UTR is 
    defined as the last 3'UTR fragment for each transcript and it will be cut 
    if there is any overlaps with other exons.
}
\usage{
    utr3Annotation(txdb, orgDbSYMBOL, MAX_EXONS_GAP = 10000)
}
\arguments{
    \item{txdb}{an object of \link[GenomicFeatures]{TxDb}}
    \item{orgDbSYMBOL}{a string indicates org SYMBOL to entriz id map}
    \item{MAX_EXONS_GAP}{maximul exon gap for distal CP site}
}
\value{
    return an object of GRanges
}
\author{
    Jianhong Ou
}
\examples{
    if(interactive()){
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
        library(org.Mm.eg.db)
        
        utr3Annotation(TxDb.Mmusculus.UCSC.mm10.knownGene, "org.Mm.egSYMBOL")
    }
}
\keyword{ misc }
