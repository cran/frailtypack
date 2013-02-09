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

\method{plot}{additivePenal}(x, ...)

}

\arguments{
  \item{x}{ An fitted additive frailty model (output from calling \code{additivePenal}}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.frailtyPenal}}}

}

\value{
     Print a plot of HR and survival function of a class \code{additivePenal} object
}

\seealso{\code{\link{print.additivePenal}},\code{\link{summary.additivePenal}},\code{\link{additivePenal}},
}


\examples{

\dontrun{

data(dataAdditive)

modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
                correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=862,
                hazard="Splines")

# It takes around 4 minutes to converge
# 'var1' is boolean as a treatment variable

plot(modAdd)

}

}

\keyword{methods}
