\name{valley}
\alias{valley}
\title{
    get the local minimal square standard error (SSE)
}
\description{
    For a giving numeric vectors, calculate the top N local minimal square standard 
    error. It will also include the 
    saved ID if it is in the range of (ss, se) 
}
\usage{
valley(x, ss, se, n = 1, savedID = NA, filterByPval = TRUE)
}
\arguments{
  \item{x}{
    numeric vector
}
  \item{ss}{
    start searching position
}
  \item{se}{
    end searching position
}
  \item{n}{
    the length of output. If n=-1, output all the local minimal SSE positions.
}
  \item{savedID}{
    saved positions
}
  \item{filterByPval}{
    logical. Filter the positions by p value or not.
}
}
\value{
    a numeric vector, position list.
}
\author{
    Jianhong Ou
}
\keyword{ misc }
