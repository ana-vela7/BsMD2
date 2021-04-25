\name{BM93e3}
\alias{BM93e3}
\docType{data}
\title{
Example 3 data in Box and Meyer (1993)
}
\description{
2^{8-4} Fractional factorial design in the injection molding example from Box, Hunter and Hunter (1978).
}
\usage{data("BM93e3")}
\format{
  A data frame with 20 observations on the following 10 variables.
  \describe{
    \item{\code{blk}}{a numeric vector. }
    \item{\code{A}}{a numeric vector. Mold temperature factor.}
    \item{\code{B}}{a numeric vector. Moisture content factor.}
    \item{\code{C}}{a numeric vector. Holding Pressure factor.}
    \item{\code{D}}{a numeric vector. Cavity thickness factor.}
    \item{\code{E}}{a numeric vector. Booster pressure factor.}
    \item{\code{F}}{a numeric vector. Cycle time factor.}
    \item{\code{G}}{a numeric vector. Gate size factor.}
    \item{\code{H}}{a numeric vector. Screw speed factor.}
    \item{\code{y}}{a numeric vector. Shrinkage response.}
  }
}
\source{
Box G. E. P, Hunter, W. C. and Hunter, J. S. (1978). Statistics for Experimenters. Wiley.

Box G. E. P, Hunter, W. C. and Hunter, J. S. (2004). Statistics for Experimenters II. Wiley.

Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors in Fractionated Screening Experiments". Journal of Quality Technology. Vol. 25. No. 2. pp. 94–105.
}
\examples{
library(BsMD2)
data(BM93e3)
print(BM93e3)
}
\keyword{datasets}
