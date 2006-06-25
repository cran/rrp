\name{XPtrToDist}
\alias{XPtrToDist}
\title{Coercion from an XPtr object into a dist one}
\description{
Coercies  an \code{XPtr} object into a \code{dist} object which can be
used in other algorithms.
}
\usage{
XPtrToDist(d)
}
\arguments{
  \item{d}{an \code{XPtr} object}
}
\details{
This routine does not duplicate but just makes available as a \code{dist}
 the  data contained in \code{d}
}
\value{returns an invisible object of class \code{dist}
}
\references{Iacus, S.M., Porro, G. (2006) Random Recursive Partitioning and 
its applications to missing data imputation, classification and 
average treatment effect estimation, \emph{submitted}.}
\author{S.M. Iacus}
\seealso{\code{\link{newXPtr}}}
\examples{
a <- newXPtr(10,5)
(XPtrToDist(a))
as.dist(matrix(5,10,10))
}
\keyword{tree}
