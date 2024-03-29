% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.nestedPenal.R
\name{summary.nestedPenal}
\alias{summary.nestedPenal}
\alias{print.summary.nestedPenal}
\title{summary of regression coefficient estimates of a nested frailty model}
\usage{
\method{summary}{nestedPenal}(object, level = 0.95, len = 6, d = 2,
lab="hr", ...)
}
\arguments{
\item{object}{output from a call to nestedPenal.}

\item{level}{significance level of confidence interval. Default is 95\%.}

\item{len}{the total field width. Default is 6.}

\item{d}{the desired number of digits after the decimal point. Default of 6
digits is used.}

\item{lab}{label of printed results.}

\item{\dots}{other unused arguments.}
}
\value{
Prints HR and its confidence intervals for each regression
coefficient. Confidence level is allowed (level argument).
}
\description{
This function returns hazard rations (HR) and its confidence intervals for
each regression coefficient.
}
\examples{


\donttest{

data(dataNested)

modNested <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,data=dataNested,
n.knots=8,kappa=c(50000,50000),hazard="Splines")

#- It takes 90 minutes to converge (depends on processor)

summary(modNested)

}


}
\seealso{
\code{\link{frailtyPenal}}
}
\keyword{methods}
