\name{summary.frailtyPenal}
\alias{summary.frailtyPenal}
\alias{print.summary.frailtyPenal}
\alias{intervals.frailtyPenal}

\title{summary of 'frailtyPenal'}

\description{
  This function returns hazard rations (HR) and its confidence intervals
}
\usage{
summary.frailtyPenal(object, level = 0.95, len = 6, d = 2, lab="hr", ...)
}

\arguments{
  \item{object}{output from a call to frailtyPenal.} 
  \item{level}{significance level of confidence interval. Default is 95\%.}
  \item{len}{the desired number of digits after the decimal point. Default of 6 digits is used.}
  \item{d}{the total field width. Default is 6.}
  \item{lab}{label of printed results.}
  \item{\dots}{ other unusued arguments.}
}

\details{
  This function calls to intervals.frailtyPenal 
}

\value{
  Prints HR and its confidence intervals. Confidence level is allowed (level argument) 
}


\seealso{ \code{\link{intervals.frailtyPenal}}}
\examples{

data(kidney)
mod<-frailtyPenal(Surv(time,status)~age+sex+cluster(id),
     data=kidney,n.knots=8,kappa1=1000)
summary(mod)

# confidence interval at 99% level 

summary(mod,level=0.99)

}


\keyword{ survival }

