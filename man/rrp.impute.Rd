\name{rrp.impute}
\alias{rrp.impute}
\title{Nearest neighbor hot-deck imputation using RRP dissimilarity matrix}
\description{
This function performs a simple nearest neighbor hot-deck imputation method
using  the RRP dissimilarity matrix.
}
\usage{
rrp.impute(data, D = NULL, k = 1, msplit = 10, Rep = 250, cut.in = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a \code{data.frame} containing missing data on some covariates}
  \item{D}{\code{NULL} or a dist object with attribute method = `RRP'}
  \item{k}{number of nearest neighbors to use}
  \item{msplit}{minimum split parameter in the \code{rpart} algorithm}
  \item{Rep}{number of RRP replications}
  \item{cut.in}{number of breaks used to cut continuous covariates}
}
\details{
If missing data are on a continuous covariate, the missing value is imputed
as the average of the covariate values of the nearest neighbors, otherwise
the majority of the `votes' determines the class of the missing observation
on the basis of nearest available data.

If \code{D} is \code{NULL} a RRP-dissimilarity matrix is created.
}
\value{
A list
  \item{new.data}{a copy of the data \code{data} with missing data imputed}
  \item{dist}{an object of class \code{dist} which is a copy of the RRP 
  dissimilarity matrix generated or passed in input and used to search 
  for nearest neighbors}
}
\references{
Iacus, S.M., Porro, G. (2006) Random Recursive Partitioning and 
its applications to missing data imputation, classification and 
average treatment effect estimation, \emph{submitted}.
}
\author{S.M. Iacus}
\seealso{\code{\link{rrp.dist}}, \code{\link{rrp.class}}}
\examples{
data(iris)

X <- iris
n <- dim(X)[1]

set.seed(123)
miss <- sample(1:n, 10)
for(i in miss)
 X[i, sample(1:5, 2)] <- NA
 
X[miss,] 

## unsupervised
x <- rrp.impute(X)

x$new.data[miss,]
iris[miss,]
}
\keyword{multivariate}
\keyword{datagen}
