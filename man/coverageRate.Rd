\name{coverageRate}
\alias{coverageRate}
\title{
    coverage rate of genes and 3UTRs
}
\description{
    calculate coverage rate of gene and 3UTRs. This function is used for quality control of the coverage.
    The coverage rate can show the complexity of RNA-seq library.
}
\usage{
    coverageRate(coverage, txdb, genome, 
                         cutoff_readsNum=1,
                         cutoff_expdGene_cvgRate=0.1,
                         cutoff_expdGene_sampleRate=0.5,
                         which=NULL, ...)
}
\arguments{
    \item{coverage}{coverage for each sample, 
        output of \link{coverageFromBedGraph}}
    \item{txdb}{an object of \link[GenomicFeatures]{TxDb}}
    \item{genome}{an object of BSgenome}
    \item{cutoff_readsNum}{cutoff reads number. If the coverage in the location 
    is greater than cutoff_readsNum, the location will be treated as coveraged 
    by signal.}
    \item{cutoff_expdGene_cvgRate,cutoff_expdGene_sampleRate}{
    cutoff_expdGene_cvgRate and cutoff_expdGene_sampleRate are the parameters used
    to calculate which gene is expressed in all input dataset. 
    cutoff_expdGene_cvgRate set the cutoff value for the coverage rate of each gene;
    cutoff_expdGene_sampleRate set the cutoff value for ratio of numbers of expressed 
    and all samples for each gene.
    for example, by default, cutoff_expdGene_cvgRate=0.1 and cutoff_expdGene_sampleRate=0.5,
    surppose there are 4 samples, for one gene, 
    if the coverage rates by base are: 0.05, 0.12, 0.2, 0.17, 
    this gene will be count as expressed gene because 
        mean(c(0.05, 0.12, 0.2, 0.17) > cutoff_expdGene_cvgRate) > cutoff_expdGene_sampleRate
    if the coverage rates by base are: 0.05, 0.12, 0.07, 0.17,
    this gene will be count as un-expressed gene because
        mean(c(0.05, 0.12, 0.07, 0.17) > cutoff_expdGene_cvgRate) <= cutoff_expdGene_sampleRate
    }
    \item{which}{
        an object of \link[GenomicRanges]{GRanges} or NULL. If it is not NULL,
        only the exons overlapping the given ranges are used.
    }
    \item{...}{not used.}
}
\value{
    return a datafrom with colnames : 
        gene.coverage.rate: coverage per base for all genes, 
        expressed.gene.coverage.rate: coverage per base for expressed genes, 
        UTR3.coverage.rate: coverage per base for all 3' UTRs,
        UTR3.expressed.gene.subset.coverage.rate: coverage per base for 3' UTRs of expressed genes.
    and rownames:
        the names of coverage.
}
\author{
    Jianhong Ou
}
\examples{
    if(interactive()){
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        path <- file.path(find.package("InPAS"), "extdata")
        bedgraphs <- c(file.path(path, "Baf3.extract.bedgraph"), 
                        file.path(path, "UM15.extract.bedgraph"))
        hugeData <- FALSE
        
        coverage <- coverageFromBedGraph(bedgraphs, 
                                 tags=c("Baf3", "UM15"),
                                 genome=BSgenome.Mmusculus.UCSC.mm10,
                                 hugeData=hugeData)
        coverageRate(coverage,
                    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,  
                    genome=BSgenome.Mmusculus.UCSC.mm10,
                    which = GRanges("chr6", ranges=IRanges(98013000, 140678000)))
    }
}
\keyword{ misc }
