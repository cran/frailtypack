
\name{plot.epoce}
\alias{plot.epoce}

\title{Plot values of estimators of the Expected Prognostic Observed Cross-Entropy (EPOCE).}

\description{
  Plots values of estimators MPL and CVPL for evaluating EPOCE. No confidence interval.
}

\usage{
  \method{plot}{epoce}(x, pos.legend="topright", cex.legend=0.7, ...)
}

\arguments{
  \item{x}{An object from the 'prediction' function, more generaly a \code{predFrailty} class object.}
  \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list
                    '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and '"center"'. The default is '"topright"'.}
  \item{cex.legend}{size of the legend. Default is 0.7.}
  \item{\dots}{Other unused arguments.}
}

\value{
  Print one plot with two curves. One for the Mean Prognosis Loss (MPL) and an other for the Cross-Validated Prognosis Loss (CVPL).
}

\keyword{file}
