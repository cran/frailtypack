% --- Source file: plot.nestedPenal.Rd ---
\name{plot.nestedPenal}
\Rdversion{1.1}
\alias{plot.nestedPenal}
\alias{lines.nestedPenal}


\title{Plot Method for a Nested frailty model.}

\description{
  Plots estimated baseline survival and hazard functions (output from an object of class `frailtyPenal' for nested frailty models). Confidence bands are allowed. 
}

\usage{

\method{plot}{nestedPenal}(x, type.plot = "hazard", conf.bands=TRUE, ...)

}

\arguments{
  \item{x}{ A nested model, more generaly an object of class \code{frailtyPenal for Nested frailty models} (output from calling \code{frailtyPenal}
    function).}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted.  The
                   default is to do so.}
  \item{\dots}{ Other graphical parameters}
}

\value{
      Print a plot of the baseline survival or hazard functions with the confidence bands or not (conf.bands argument)
}

\seealso{\code{\link{print.nestedPenal}},\code{\link{summary.nestedPenal}},\code{\link{frailtyPenal for Nested frailty models}}}


\examples{

\dontrun{

  data(dataNested)
  modNested<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+subcluster(subgroup)+
  cov1+cov2,data=dataNested,n.knots=8,kappa1=50000)

# It takes around 90 minutes to converge (depends on the processor)  #

  plot(modNested,conf.bands=FALSE)

}


}

\keyword{methods}
