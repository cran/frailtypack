\name{plot.frailtyPenal}
\alias{plot.frailtyPenal}
\alias{lines.frailtyPenal}


\title{Plot Method for an object of class `frailtyPenal'.}

\description{
  Plots estimated baseline survival and hazard functions from an object of class `frailtyPenal'. Confidence 
  bands are allowed. 
}

\usage{
plot.frailtyPenal(x, type.plot = "hazard", conf.bands=TRUE, ...)
}

\arguments{
  \item{x}{ Object of class \code{frailtyPenal} (output from calling \code{frailtyPenal}
    function).}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted.  The
                   default is to do so.}
  \item{\dots}{ Other graphical parameters}
}

\value{
     Print a plot of class \code{frailtyPenal}
}

\seealso{\code{\link{print.frailtyPenal}}}


\examples{

data(readmission)

# Let's compare shared frailty model with Cox proportional hazards model 
mod.sha<-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
             n.knots=10,kappa1=10000,data=readmission)
plot(mod.sha,type="surv",conf=FALSE)
mod.cox<-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
             n.knots=10,kappa1=10000,data=readmission,Frailty=FALSE)
lines(mod.cox,type="surv",conf=FALSE,col=2)




# Stratified model
mod<-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
             n.knots=10,kappa1=10000,kappa2=10000,data=readmission)
plot(mod)

# no confidence bands
plot(mod,conf.bands=FALSE)

}

\keyword{ survival }
