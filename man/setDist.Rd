\name{setDist}
\alias{setDist}
\title{Sets the elements of a dist object to a constant vector-wise}
\description{
Sets the elements of a dist object to a constant vector-wise. The idea is that 
a symmetrix matrix \code{M} can be represented as a triangular object {d}. 
With this function you set \code{M[x,x] = k} but instead of working on 
\code{M} this function works on its \code{dist} representation \code{d}.
}
\usage{
setDist(d, x, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{a dist object, upper diagonal, representing a symmetric matrix}
  \item{x}{a vector of indices}
  \item{k}{the constant to the elements will be set to}
}
\value{
  returns the modified dist object
}
\author{S. M. Iacus}
\seealso{\code{\link{dist}}}
\examples{
M <- matrix(1:100,5,20)
M <- M %*% t(M)
d <- as.dist(M)
x <- 1:3

setDist(d,x,0)
d
# which is the equivalent of (apart for the diagonal elements)
M[x,x] <- 0

# Indeed, we have
d
as.dist(M)
}
\keyword{misc}
