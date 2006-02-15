\name{rrp.dist}
\alias{rrp.dist}
\title{RRP dissimilarity matrix}
\description{
Main piece of the RRP algorithm which generates the RRP dissimilarity
matrix.
}
\usage{
rrp.dist(X, treated = NULL, msplit = 10, Rep = 250, cut.in = 15, 
 check.bal = FALSE, plot = FALSE)
}
\arguments{
  \item{X}{a \code{data.frame} object}
  \item{treated}{optional class indicator variable}
  \item{msplit}{minimum split parameter in the \code{rpart} algorithm}
  \item{Rep}{number of RRP replications}
  \item{cut.in}{number of breaks in which to cut continuous variables}
  \item{check.bal}{indicator function. If \code{TRUE} balance check using
  hyper-rectangles will be used inside leaves}
  \item{plot}{wheter to plot the porximity matrix as image}
}
\details{
This algorithm allows for missing data in \code{X}.
}
\value{an object of class \code{dist}
}
\references{
Iacus, S.M., Porro, G. (2006) Random recursive partitioning: a 
matching method for the estimation of the average treatment effect
}
\author{S.M. Iacus}
\seealso{\code{\link{rpart}}}
\keyword{tree}
