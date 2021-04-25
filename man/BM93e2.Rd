\name{BM93e2}
\alias{BM93e2}
\docType{data}
\title{
Example 2 data in Box and Meyer (1993)
}
\description{
12-run Plackett-Burman design for the study of fatigue life of weld repaired castings.
}
\usage{data("BM93e2")}
\format{
  A data frame with 12 observations on the following 8 variables.
  \describe{
    \item{\code{A}}{a numeric vector. Initial structure factor.}
    \item{\code{B}}{a numeric vector. Bead size factor.}
    \item{\code{C}}{a numeric vector. Pressure treat factor.}
    \item{\code{D}}{a numeric vector. Heat treat factor.}
    \item{\code{E}}{a numeric vector. Cooling rate factor.}
    \item{\code{F}}{a numeric vector. Polish factor.}
    \item{\code{G}}{a numeric vector. Final treat factor.}
    \item{\code{y}}{a numeric vector. Natural log of fatigue life response.}
  }
}
\source{
Hunter, G. B., Hodi, F. S., and Eager, T. W. (1982). "High-Cycle Fatigue of Weld Repaired Cast Ti-6A1-4V". Metallurgical Transactions 13A, pp. 1589–1594.

Box, G.E.P. and R.D. Meyer (1993). "Finding the Active Factors in Fractionated Screening Experiments". Journal of Quality Technology. Vol. 25. No. 2. pp. 94-105.
}
\examples{
library(BsMD2)
data(BM93e2)
print(BM93e2)
}
\keyword{datasets}
