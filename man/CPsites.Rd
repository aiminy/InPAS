\name{CPsites}
\alias{CPsites}
\title{
    predict the cleavage and polyadenylation(CP) site
}
\description{
    predict the alternative cleavage and polyadenylation (CP or APA) site. 
}
\usage{
    CPsites(coverage, groupList=NULL, genome, utr3, 
        window_size=100, search_point_START=50, search_point_END=NA, 
        cutStart=window_size, cutEnd=0, adjust_distal_polyA_end=TRUE, 
        coverage_threshold=5, long_coverage_threshold=2, 
        background=c("same_as_long_coverage_threshold", 
                     "1K", "5K", "10K", "50K"),
        txdb=NA,
        gcCompensation=NA, mappabilityCompensation=NA, 
        FFT=FALSE, fft.sm.power=20, 
        PolyA_PWM=NA, classifier=NA, classifier_cutoff=.8, step=1, 
        shift_range=window_size,
        BPPARAM=NULL, tmpfolder=NULL, silence=TRUE)
}
\arguments{
    \item{coverage}{coverage for each sample, 
        output of \link{coverageFromBedGraph}}
    \item{groupList}{group list of tag names}
    \item{genome}{an object of \link[BSgenome]{BSgenome}}
    \item{utr3}{output of \link{utr3Annotation}}
    \item{window_size}{window size for noval distal position searching and 
        adjusted polyA searching, default: 100}
    \item{search_point_START}{start point for searching}
    \item{search_point_END}{end point for searching}
    \item{cutStart}{how many nucleotides should be removed from the start 
        before search, 0.1 means 10 percent, 25 means cut first 25.}
    \item{cutEnd}{how many nucleotides should be removed  from the end before 
        search, 0.1 means 10 percent.}
    \item{adjust_distal_polyA_end}{If true, adjust distal polyA end by 
        \link[cleanUpdTSeq]{cleanUpdTSeq}}
    \item{coverage_threshold}{cutoff coverage threshold for first 100 
        nucleotides. If the coverage of first 100 nucleotides is lower than 
        coverage_threshold, that transcript will be dropped.}
    \item{long_coverage_threshold}{cutoff threshold for coverage in the region
        of long form. If the coverage in the region of long form is less than 
        long_coverage_threshold, that transcript will be dropped.}
    \item{background}{the range for calculating cutoff threshold of local background}
    \item{txdb}{an object of \link[GenomicFeatures]{TxDb}}
    \item{gcCompensation}{GC content compensation vector}
    \item{mappabilityCompensation}{mappability compensation vector}
    \item{FFT}{use Fast Fourier Transform Algorithm to smooth the data or not. 
        default: FALSE}
    \item{fft.sm.power}{if FFT is TRUE, the frequency should be removed}
    \item{PolyA_PWM}{Position Weight Matrix of polyA}
    \item{classifier}{An object of class 
        \code{"\link[cleanUpdTSeq]{PASclassifier}"}}
    \item{classifier_cutoff}{This is the cutoff used to assign whether a 
        putative pA is true or false. 
        This can be any floating point number between 0 and 1. For example,
        classifier_cutoff = 0.5 will assign an 
        putative pA site with prob.1 > 0.5 to the True class (1), 
        and any putative pA site with prob.1 <= 0.5 as False (0).}
    \item{step}{adjust step, default 1, means adjust by each base by cleanUpdTSeq.}
    \item{shift_range}{the shift range for polyA site searching}
    \item{BPPARAM}{An optional \code{\link[BiocParallel]{BiocParallelParam}} 
        instance determining the parallel back-end to be used during 
        evaluation, or a list of BiocParallelParam instances, to be applied 
        in sequence for nested calls to bplapply.}
    \item{tmpfolder}{temp folder could save and reload the analysis data for 
    resume analysis.}
    \item{silence}{report progress or not. default not report.}
}
\value{
    return an object of GRanges contain the estimated CP sites.
}
\references{
    ref: Cheung MS, Down TA, Latorre I, Ahringer J. Systematic bias in 
        high-throughput sequencing data and its correction by BEADS. 
        Nucleic Acids Res. 2011 Aug;39(15):e103. doi: 10.1093/nar/gkr425. 
        Epub 2011 Jun 6. PubMed PMID: 21646344;
        PubMed Central PMCID: PMC3159482.

    mappability could be calculated by 
        [GEM](http://algorithms.cnag.cat/wiki/Man:gem-mappability)
 
    ref: Derrien T, Estelle J, Marco Sola S, Knowles DG, Raineri E, Guigo R, 
        Ribeca P. 
        Fast computation and applications of genome mappability. PLoS One.
        2012;7(1):e30377. doi: 10.1371/journal.pone.0030377. 
        Epub 2012 Jan 19. PubMed
        PMID: 22276185; PubMed Central PMCID: PMC3261895.
}
\author{
    Jianhong Ou
}
\examples{
    if(interactive()){
        library(BSgenome.Mmusculus.UCSC.mm10)
        path <- file.path(find.package("InPAS"), "extdata")
        bedgraphs <- file.path(path, "Baf3.extract.bedgraph")
        data(utr3.mm10)
        tags <- "Baf3"
        genome <- BSgenome.Mmusculus.UCSC.mm10
        coverage <- 
            coverageFromBedGraph(bedgraphs, tags, genome, hugeData=FALSE)
        CP <- CPsites(coverage=coverage, gp1=tags, gp2=NULL, genome=genome, 
            utr3=utr3.mm10, coverage_threshold=5, long_coverage_threshold=5)
    }
}
\keyword{ misc }
