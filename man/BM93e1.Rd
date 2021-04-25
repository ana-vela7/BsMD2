\name{BM93e1}
\alias{BM93e1}
\docType{data}
\title{
Example 1 data in Box and Meyer (1993)
}
\description{
12-run Plackett-Burman design from the 2^5 reactor example from Box, Hunter and Hunter (1977).
}
\usage{data("BM93e1")}
\format{
  A data frame with 12 observations on the following 7 variables.
  \describe{
    \item{\code{Run}}{a numeric vector. Run number from a 2^5 factorial design in standard order.}
    \item{\code{A}}{a numeric vector. Feed rate factor.}
    \item{\code{B}}{a numeric vector. Catalyst factor.}
    \item{\code{C}}{a numeric vector. Agitation factor.}
    \item{\code{D}}{a numeric vector. Temperature factor.}
    \item{\code{E}}{a numeric vector. Concentration factor.}
    \item{\code{y}}{a numeric vector. Percent reacted response.}
  }
}
\source{
Box, G.E.P, Hunter, W.C. and Hunter, J. S. (1978). Statistics for Experimenters. Wiley.

Box, G.E.P. and R.D. Meyer (1993). "Finding the Active Factors in Fractionated Screening Experiments". Journal of Quality Technology. Vol. 25. No. 2. pp. 94-105.
}
\examples{
library(BsMD2)
data(BM93e1)
print(BM93e1)
}
\keyword{datasets}
