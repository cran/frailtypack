
% --- Source file: plot.jointPenal.Rd ---
\name{plot.jointPenal}
\Rdversion{1.1}
\alias{plot.jointPenal}
\alias{lines.jointPenal}


\title{Plot Method for a Joint frailty model.}

\description{
  Plots estimated baseline survival and hazard functions of a joint frailty model (output from an object of class 'frailtyPenal' for joint frailty models ) for each type of event (terminal or recurrent). Confidence bands are allowed. 
}

\usage{

\method{plot}{jointPenal}(x, event = "both", type.plot = "hazard", conf.bands
                 = FALSE, ...)


}

\arguments{
  \item{x}{ A joint model, or more generaly an object of class \code{frailyPenal for Joint frailty model} (output from calling \code{frailtyPenal}
    function).}
 \item{event}{ a character string specifying the type of curve. Possible
                   value are "terminal", "recurrent", or "both". The default
                   is "both".}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{conf.bands}{ logical value. Determines whether confidence bands will be plotted. The default is to do so.}
  \item{\dots}{ Other graphical parameters}
}

\value{
     Print a plot of the baseline survival or hazard functions for each type of event or both with the confidence bands or not (conf.bands argument)
}

\seealso{\code{\link{print.jointPenal}},\code{\link{summary.jointPenal}},\code{\link{frailtyPenal for Joint frailty models}}

}


\examples{


\dontrun{

  data(dataJoint)

  modJoint<-frailtyPenal(Surv(time.entry,time.end,status)~cluster(id)+var1+var2
                         +terminal(status.terminal),
                         formula.terminalEvent=~var1,
                         data=dataJoint,n.knots=7,
                         kappa1=1, kappa2=1,Frailty=TRUE, joint=TRUE, recurrentAG=TRUE)

# It takes around 1 minute to converge. #


plot(modJoint,type.plot="haz",event="recurrent",conf.bands=TRUE)
plot(modJoint,type.plot="haz",event="terminal",conf.bands=TRUE)
plot(modJoint,type.plot="haz",event="both",conf.bands=TRUE))
plot(modJoint,type.plot="su",event="recurrent",conf.bands=TRUE))
plot(modJoint,type.plot="su",event="terminal",conf.bands=TRUE))
plot(modJoint, event="terminal",type.plot="su",conf.bands=TRUE))
plot(modJoint,type.plot="su","terminal",conf.bands=TRUE))

}


}

\keyword{methods}
