% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.jointPenal.R
\name{plot.jointPenal}
\alias{plot.jointPenal}
\alias{lines.jointPenal}
\title{Plot Method for a Joint frailty model.}
\usage{
\method{plot}{jointPenal}(x, event = "Both", type.plot = "Hazard", conf.bands
= FALSE, pos.legend="topright", cex.legend = 0.7, ylim, main, color = 2, median=TRUE,
Xlab = "Time", Ylab = "Hazard function", ...)
}
\arguments{
\item{x}{A joint model, i.e. an object of class \code{frailtyPenal} for
Joint frailty model (output from calling \code{frailtyPenal} function).}

\item{event}{a character string specifying the type of curve. Possible value
are "Terminal", "Recurrent", or "Both". The default is "Both".}

\item{type.plot}{a character string specifying the type of curve. Possible
value are "Hazard", or "Survival". The default is "Hazard". Only the first
letters are required, e.g "Haz", "Su"}

\item{conf.bands}{logical value. Determines whether confidence bands will be
plotted. The default is to do so.}

\item{pos.legend}{The location of the legend can be specified by setting
this argument to a single keyword from the list '"bottomright"', '"bottom"',
'"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
'"center"'. The default is '"topright"'}

\item{cex.legend}{character expansion factor *relative* to current
'par("cex")'. Default is 0.7}

\item{ylim}{y-axis limits}

\item{main}{plot title}

\item{color}{curve color (integer)}

\item{median}{Logical value. Determines whether survival median will be plotted. Default is TRUE.}

\item{Xlab}{Label of x-axis. Default is '"Time"'}

\item{Ylab}{Label of y-axis. Default is '"Hazard function"'}

\item{...}{other unused arguments}
}
\value{
Print a plot of the baseline survival or hazard functions for each
type of event or both with the confidence bands or not (conf.bands argument)
}
\description{
Plots estimated baseline survival and hazard functions of a joint frailty
model (output from an object of class 'JointPenal' for joint frailty models
) for each type of event (terminal or recurrent). Confidence bands are
allowed.
}
\examples{


\donttest{

data(readmission)

#-- Gap-time
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+
charlson+terminal(death),formula.terminalEvent=~sex+dukes+charlson,
data=readmission,n.knots=14,kappa=c(100,100))

#-- It takes around 1 minute to converge --#

plot(modJoint.gap,type.plot="Haz",event="recurrent",conf.bands=TRUE)
plot(modJoint.gap,type.plot="Haz",event="terminal",conf.bands=TRUE)
plot(modJoint.gap,type.plot="Haz",event="both",conf.bands=TRUE)

plot(modJoint.gap,type.plot="Su",event="recurrent",conf.bands=TRUE)
plot(modJoint.gap,type.plot="Su",event="terminal",conf.bands=TRUE)
plot(modJoint.gap,type.plot="Su",event="both",conf.bands=TRUE)


}


}
\seealso{
\code{\link{frailtyPenal}}
}
\keyword{methods}
