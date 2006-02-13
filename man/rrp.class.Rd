\name{rrp.class}
\alias{rrp.class}
\title{Nonparametric nearest neighbor classification unsing RRP dissimilarity}
\description{
This method is a simple nearest neighbor based on the RRP dissimilarity matrix.
It can be used also as a supervised method, i.e. using the class variable
in the construction of the RRP dissimilarity matrix and by imputing missing 
values to the test set.
}
\usage{
rrp.class(x, cl, train, test, k = 1)
}
\arguments{
  \item{x}{a \code{dist} object with attribute \code{method} = \code{RRP}}
  \item{cl}{class vector, coerced to be of type \code{factor}}
  \item{train}{the vecotr of indexes of the training set}
  \item{test}{the vector of training indexes of test set}
  \item{k}{number of nearest to consider}
}
\value{
a vector of type \code{factor} with predicted classes.
}
\references{Iacus, S.M., Porro, G. (2006) Random Recursive Partitioning and 
its applications to missing data imputation, classification and 
average treatment effect estimation, \emph{submitted}.}
\author{S.M. Iacus}
\seealso{\code{\link{rrp.dist}}}
\examples{
data(iris)
X <- iris[,-5]
n <- dim(X)[1]

set.seed(123)
test <- sample(1:n, 10)
train <- (1:n)[-test]

## unsupervised
D <- rrp.dist(X)
pred <- rrp.class(D, iris[train,5], train, test)
table(pred, iris[test, 5])

# supervised
X <- iris
X[test,5] <- NA
D <- rrp.dist(X)
pred <- rrp.class(D, iris[train,5], train, test)
table(pred, iris[test, 5])
}
\keyword{classif}
