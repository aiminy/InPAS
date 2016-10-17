\name{trimSeqnames}
\alias{trimSeqnames}
\title{
    trim the sequence names
}
\description{
    only ^chr[0-9XY]+$ is OK.
}
\usage{
trimSeqnames(genome)
}
\arguments{
  \item{genome}{
    an BSgenome object
}
}
\value{
an character vector with trimmed seqnames
}
\author{
    Jianhong Ou
}
\keyword{ misc }