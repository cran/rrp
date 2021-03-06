\name{DWvsPSID}
\alias{DWvsPSID}
\docType{data}
\title{Dehejia-Wahba treated units versus PSID control individuals}
\description{
The Lalonde subset of DW treated units versus PSID (Panel Study of Income 
Dynamics) control individuals
}
\usage{data(DWvsPSID)}
\format{
  A data frame with 2675 observations on the following 10 variables.
  \describe{
    \item{\code{treated}}{treated variable indicator}
    \item{\code{age}}{age}
    \item{\code{education}}{years of education}
    \item{\code{black}}{race indicator variable}
    \item{\code{married}}{marital status indicator variable}
    \item{\code{nodegree}}{indicator variable of not possessing a degree}
    \item{\code{re74}}{real earnings in 1974}
    \item{\code{re75}}{real earnings in 1975}
    \item{\code{re78}}{real earnings in 1978 (post treatment outcome)}
    \item{\code{hispanic}}{ethnic indicator variable}
  }
}
\details{
These two sets of treated and control units can be hardly matched.
}
\source{see references}
\references{
Dehejia, R., Wahba, S. (1999) Causal Effects in Nonexperimental Studies: 
Reevaluating the Evaluation of Training Programs,  \emph{Journal of 
the American Statistical Association}, 94, 1053-1062.
}
\keyword{datasets}
