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

\method{plot}{additivePenal}(x, type.plot = "hazard", conf.bands=TRUE, ...)

}

\arguments{
  \item{x}{ An fitted additive frailty model (output from calling \code{additivePenal}
    function).}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted.  The
                   default is to do so.}
  \item{\dots}{ Other graphical parameters}
}

\value{
     Print a plot of HR and survival function of a class \code{additivePenal} object
}

\seealso{\code{\link{print.additivePenal}},\code{\link{summary.additivePenal}},\code{\link{additivePenal}},
}


\examples{


\dontrun{
  
  data(dataAdditive) 
  modAdd<-additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
                 correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=862)

# It takes around 4 minutes to converge. 'var1' is boolean as a treatment variable. #

  plot(modAdd)

} 



}

\keyword{methods}
