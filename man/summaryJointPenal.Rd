% --- Source file: summaryjointPenal.Rd ---
\name{summary.jointPenal}
\Rdversion{1.1}
\alias{summary.jointPenal}
\alias{print.summary.jointPenal}


\title{summary of parameter estimates of a joint frailty model}

\description{
  This function returns hazard rations (HR) and its confidence intervals
}
\usage{
\method{summary}{jointPenal}(object, level = 0.95, len = 6, d = 2, lab="hr", ...)
}

\arguments{
  \item{object}{output from a call to frailtyPenal for joint models} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unusued arguments.}
}


\value{
  Prints HR and its confidence intervals for each covariate. Confidence level is allowed (level argument) 
}


\seealso{ \code{\link{frailtyPenal for Joint frailty models}},
\code{\link{plot.jointPenal}},
\code{\link{print.jointPenal}},
}
\examples{

\dontrun{

  data(readmission)

## Gap-time ##
modJoint_gap<-frailtyPenal(Surv(time,event)~cluster(id)+sex+as.factor(dukes)
	      +as.factor(charlson)+terminal(death),
	      formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
	      data=readmission,n.knots=14,kappa1=9550000000,kappa2=1410000000000,
	      Frailty=TRUE,joint=TRUE,recurrentAG=FALSE,hazard="Splines")

## Calendar time ##
modJoint_calendar<-frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
	      +as.factor(dukes)+as.factor(charlson)+terminal(death),
	      formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
	      data=readmission,n.knots=10,kappa1=9550000000,
	      kappa2=1410000000000,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE
	      ,hazard="Splines")
	 
# It takes around 1 minute to converge  

  summary(modJoint_gap)
  summary(modJoint_calendar)
}


}


\keyword{methods}
