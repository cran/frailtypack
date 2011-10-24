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

\method{plot}{nestedPenal}(x, ...)

}

\arguments{
  \item{x}{ A nested model, more generaly an object of class \code{frailtyPenal for Nested frailty models} (output from calling \code{frailtyPenal}
    function).}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.frailtyPenal}}}
}

\value{
      Print a plot of the baseline survival or hazard functions with the confidence bands or not (conf.bands argument)
}

\seealso{\code{\link{print.nestedPenal}},\code{\link{summary.nestedPenal}},\code{\link{frailtyPenal for Nested frailty models}}}


\examples{

\dontrun{

  data(dataNested)
  modNested<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+
          subcluster(subgroup)+cov1+cov2,data=dataNested,
	  n.knots=8,kappa1=50000,hazard="Splines")

# It takes around 90 minutes to converge (depends on the processor)  #

  plot(modNested,conf.bands=FALSE)

}


}

\keyword{methods}
