\name{summary.multivePenal}
\Rdversion{1.1}
\alias{summary.multivePenal}

\title{
summary of parameter estimates of a multivariate frailty model.
}
\description{
  This function returns hazard ratio (HR) and its confidence intervals.
}
\usage{
\method{summary}{multivePenal}(object, level = 0.95, len = 6, d = 2, lab = "hr", ...)
}

\arguments{
  \item{object}{output from a call to multivePenal for multive joint models} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{d}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{len}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{other unusued arguments.}
}

\value{
  Prints HR and its confidence intervals for each covariate. Confidence level is allowed (level argument) 
}

\seealso{
 \code{\link{multivePenal}}
}

\examples{
\dontrun{

data(dataMultiv)
dataMultiv$timegap<-dataMultiv$t1-dataMultiv$t0

## Calendar-time ##

MultivSpli <-multivePenal(Surv(timegap,deltaEvent1)~cluster(id)+var1+var2+terminal(deltadc)+event2(deltaEvent2),formula.terminalEvent=~var1,formula2=~var1+var2+var3,data=dataMultiv,Frailty=TRUE, recurrentAG=F,cross.validation=T,n.knots=c(4,4,4),kappa=c(1,1,1),maxit=20,hazard="Splines",initialize=T)

## summary a fit
summary(MultivSpli)
}
}

\keyword{methods}
\keyword{multive}
