\name{dataJointClus}
\alias{dataJointClus}
\docType{data}
\title{
a simulated data for JointCluster.
}
\description{
a dataset.
}
\usage{data(dataJointClus)}
\format{
  A data frame with 24413 observations on the following 7 variables.
  \describe{
    \item{\code{patkey}}{an observation.}
    \item{\code{start}}{start of interval (0 or previous recurrence time)}
    \item{\code{stop}}{recurrence or censoring time}
    \item{\code{event}}{a numeric vector}
    \item{\code{cluster}}{a numeric vector}
    \item{\code{trt}}{a numeric vector}
    \item{\code{death}}{a numeric vector}
  }
}
\examples{
data(dataJointClus)
}
\keyword{datasets}
