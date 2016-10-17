## ref: Cheung MS, Down TA, Latorre I, Ahringer J. Systematic bias in 
##      high-throughput 
##      sequencing data and its correction by BEADS. Nucleic Acids Res. 2011
##      Aug;39(15):e103. doi: 10.1093/nar/gkr425. Epub 2011 Jun 6. PubMed 
##      PMID: 21646344;
##      PubMed Central PMCID: PMC3159482.

gcComp <- function(genome, seqnames, window=50){
    fref <- sapply(seqnames, function(.ele) 
        GC(s2c(as.character(genome[[.ele]]))))
    win <- floor(window/2)
    fdat <- lapply(seqnames, function(seq, ref){
        s <- s2c(as.character(genome[[seq]]))
        N <- rep("N", win)
        ss <- win+1:length(s)
        se <- ss+win
        s <- c(N, s, N, "N", "N")
        f <- mapply(function(.s, .e) GC(s[.s:.e]), ss, se)
        wk <- fref[seq]/f
    })
}

## mappability is calculated by 
##      [GEM](http://algorithms.cnag.cat/wiki/Man:gem-mappability)
## ref: Derrien T, Estellé J, Marco Sola S, Knowles DG, Raineri E, 
##      Guigó R, Ribeca P. 
##      Fast computation and applications of genome mappability. PLoS One.
##      2012;7(1):e30377. doi: 10.1371/journal.pone.0030377. Epub 2012 Jan 19. 
##      PubMed PMID: 22276185; PubMed Central PMCID: PMC3261895.
##export PATH=$PATH:~/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
##./gem-indexer -i \
## genome.fa \
## -o mm10.index.gem
##./gem-mappability -I mm10.index.gem.gem -l 100 -o mm10
##./gem-2-wig -I mm10.index.gem.gem -i mm10.mappability -o mm10.mappability.wig

mapComp <- function(mi){
    max(mi)/mi
}
