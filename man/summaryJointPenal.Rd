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

  data(dataJoint)

  modJoint<-frailtyPenal(Surv(time.entry,time.end,status)~cluster(id)+var1+var2
                         +terminal(status.terminal),
                         formula.terminalEvent=~var1,
                         data=dataJoint,n.knots=7,Frailty=TRUE,
                         kappa1=1, kappa2=1, joint=TRUE, recurrentAG=TRUE)

# It takes around 1 minute to converge  

  summary(modJoint)

}


}


\keyword{methods}
