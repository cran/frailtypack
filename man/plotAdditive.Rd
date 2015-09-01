% --- Source file: plotadditivePenal.Rd ---
\name{plot.additivePenal}
\Rdversion{1.1}
\alias{plot.additivePenal}
\alias{lines.additivePenal}


\title{Plot Method for an Additive frailty model.}

\description{
  Plots estimated baseline survival and hazard functions of an additive frailty model, more generally of a class `additivePenal' object. Confidence bands are allowed. 
}

\usage{

\method{plot}{additivePenal}(x, type.plot="hazard", conf.bands=TRUE,
 pos.legend="topright", cex.legend=0.7, main, color=2, ...)

}

\arguments{
  \item{x}{ A fitted additive frailty model (output from calling \code{additivePenal})}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted. The default is to do so.}
  \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
     '"center"'. The default is '"topright"'}
  \item{cex.legend}{character expansion factor *relative* to current
          'par("cex")'. Default is 0.7}
  \item{main}{plot title}
  \item{color}{curve color (integer)}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.frailtyPenal}}}

}

\value{
     Print a plot of HR and survival function of a class \code{additivePenal} object
}

\seealso{
 \code{\link{additivePenal}}
}


\examples{

\dontrun{

data(dataAdditive)

modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
correlation=TRUE,data=dataAdditive,n.knots=8,kappa=862,hazard="Splines")

#-- 'var1' is boolean as a treatment variable

plot(modAdd)

}

}

\keyword{methods}
