\name{utr3.hg19}
\alias{utr3.hg19}
\docType{data}
\title{
    3'UTR annotation for hg19 obtained from utr3Annotation
}
\description{
    3'UTR annotation obtained from utr3Annotation by 
    TxDb.Hsapiens.UCSC.hg19.knownGene and org.Hs.eg.db
}
\usage{data(utr3.hg19)}
\format{
    GRanges with slot start holding the start position of the 3'UTR, slot end 
    holding the end position of the 3'UTR, slot names holding transcripts and 
    gene names of 3'UTR, slot seqnames holding the chromosome location where 
    the 3'UTR is located and slot strand for strand of 3'UTR. In addition, the 
    following variables are included.
    \describe{
        \item{\code{feature}}{should be unknown or proximalCP_XXXXXXXX}
    
        \item{\code{id}}{should be utr3 or next.exon.gap}
    
        \item{\code{exon}}{exon id}
    
        \item{\code{transcript}}{trnscript id}
    
        \item{\code{gene}}{entriz gene id}
    
        \item{\code{symbol}}{gene symbol}
    }
}
\details{
    used in the examples
    Annotation data obtained by:
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    library(org.Hs.eg.db)
        
    utr3Annotation(TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.egSYMBOL)
}
\value{
    an object of GRanges.
}
\examples{
    data(utr3.hg19)
    head(utr3.hg19)
}
\keyword{datasets}
