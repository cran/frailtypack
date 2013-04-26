\name{plot.multivePenal}
\alias{plot.multivePenal}

\title{
Plot Method for a multivariate frailty model.
}
\description{
  Plots estimated baseline survival and hazard functions of a multivariate frailty model (output from an object of class 'multivePenal' for multivariate frailty models ) for each type of event (recurrent, terminal and second recurrent). Confidence bands are allowed. 
}
\usage{
\method{plot}{multivePenal}(x, event = "both", type.plot = "hazard", conf.bands = FALSE,
                            pos.legend = "topright", cex.legend = 0.7, ylim, main, ...)
}
\arguments{
  \item{x}{ A joint model, or more generaly an object of class \code{frailyPenal for Joint frailty model} (output from calling \code{frailtyPenal}
    function).}
 \item{event}{ a character string specifying the type of curve. Possible
                   value are "terminal", "recurrent", "recurrent2", or "both". The default
                   is "both".}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted. The default is to do so.}
 \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
     '"center"'. The default is '"topright"'}
  \item{cex.legend}{character expansion factor *relative* to current
          'par("cex")'. Default is 0.7}
  \item{ylim}{y-axis limits}
  \item{main}{plot title}
  \item{\dots}{ Other graphical parameters}
}


\value{
     Print a plot of the baseline survival or hazard functions for each type of event or both with the confidence bands or not (conf.bands argument)
}

\seealso{
 \code{\link{multivePenal}}
}

\examples{
\dontrun{
data(dataMultiv)

## Calendar-time ##

data(dataMultiv)
dataMultiv$timegap<-dataMultiv$t1-dataMultiv$t0

MultivSpli <- multivePenal(Surv(timegap,deltaEvent1)~cluster(id)+var1+var2+terminal(deltadc)+event2(deltaEvent2),
formula.terminalEvent=~var1,formula2=~var1+var2+var3,data=dataMultiv,Frailty=TRUE,recurrentAG=F,
cross.validation=T,n.knots=c(4,4,4),kappa=c(1,1,1),hazard="Splines",initialize=T)

## plot
plot(MultivSpli,type.plot="haz",event="recurrent",conf.bands=TRUE)
}
}
\keyword{methods}
