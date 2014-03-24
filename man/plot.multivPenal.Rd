\name{plot.multivPenal}
\alias{plot.multivPenal}

\title{
Plot Method for a multivariate frailty model.
}

\description{
  Plots of estimated baseline survival and hazard functions of a multivariate frailty model (output from an object of class 'multivPenal' for multivariate frailty models ) for each type of event (recurrent, terminal and second recurrent). Confidence intervals are allowed. 
}

\usage{
\method{plot}{multivPenal}(x, event = "both", type.plot = "hazard", conf.bands = FALSE,
 pos.legend = "topright", cex.legend = 0.7, ylim, main, color1="red", color2="blue", 
 colorEnd="green", ...)
}

\arguments{
  \item{x}{A joint multivariate model, or more generaly an object of class \code{multivPenal} (output from calling \code{multivPenal} function).}

 \item{event}{a character string specifying the type of outcome. Possible value are "terminal", "recurrent", "recurrent2", or "both". The default is "both".}

  \item{type.plot}{a character string specifying the type of curve. Possible value are "hazard", or "survival". The default is "hazard". Only the first words are required, e.g "haz", "su"}

  \item{conf.bands}{logical value. Determines whether confidence intervals will be plotted. The default is to do so.}

 \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
     '"center"'. The default is '"topright"'}

  \item{cex.legend}{character expansion factor *relative* to current
          'par("cex")'. Default is 0.7}

  \item{ylim}{y-axis limits}

  \item{main}{plot title}

  \item{color1}{curve color for recurrent event of type 1 (integer or color name in quotation marks)}

  \item{color2}{curve color for recurrent event of type 2 (integer or color name in quotation marks)}

  \item{colorEnd}{curve color for terminal event (integer or color name in quotation marks)}

  \item{\dots}{Other graphical parameters}
}

\value{
     Print a plot of the baseline survival or hazard functions for each type of event or both with the confidence intervals or not (conf.bands argument)
}

\seealso{
 \code{\link{multivPenal}}
}

\keyword{methods}
