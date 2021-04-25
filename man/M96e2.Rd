\name{M96e2}
\alias{M96e2}
\docType{data}
\title{
Reactor Experiment Data
}
\description{
Data of the Reactor Experiment from Box, Hunter and Hunter (1978).}
\usage{data("M96e2")}
\format{
  A data frame with 32 observations on the following 6 variables.
  \describe{
    \item{\code{A}}{a numeric vector. Feed rate factor.}
    \item{\code{B}}{a numeric vector. Catalyst factor.}
    \item{\code{C}}{a numeric vector. Agitation rate factor.}
    \item{\code{D}}{a numeric vector. Temperature factor.}
    \item{\code{E}}{a numeric vector. Concentration factor.}
    \item{\code{y}}{a numeric vector. Percentage reacted response.}
  }
}
\source{
Box G. E. P, Hunter, W. C. and Hunter, J. S. (2004). Statistics for Experimenters II. Wiley.

Box G. E. P, Hunter, W. C. and Hunter, J. S. (1978). Statistics for Experimenters. Wiley.
}
\examples{
library(BsMD2)
data(M96e2)
print(M96e2)
}
\keyword{datasets}
