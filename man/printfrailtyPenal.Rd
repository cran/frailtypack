\name{print.frailtyPenal}

\Rdversion{1.1}
\alias{print.frailtyPenal}

\title{Print a Short Summary of parameter estimates of a shared gamma frailty model}
\description{
  Prints a short summary of parameter estimates of a 'frailtyPenal' object 
}
\usage{
\method{print}{frailtyPenal}(x, digits = max(options()$digits - 4, 6), ...)
}
\arguments{
  \item{x}{the result of a call to the frailtyPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\value{
 

 \item{n}{the number of observations used in the fit.}
 \item{n.groups}{the maximum number of groups used in the fit}
 \item{n.events}{the number of events observed in the fit}

 \item{theta}{variance of frailty parameter}
 \item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }

  \item{SE(H)}{the standard error of the estimates deduced from the variance matrix of theta and of the coefficients.}
 \item{SE(HIH)}{the standard error of the estimates deduced from the robust estimation of the variance matrix of theta and of the coefficients.}
 \item{p}{p-value}


}


\seealso{
 \code{\link{summary.frailtyPenal}},  
\code{\link{frailtyPenal for Shared frailty models}},
 \code{\link{plot.frailtyPenal}}
  }

\examples{

data(readmission)

# /*** Shared frailty model ***/

modSha<-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
             n.knots=10,kappa1=10000,data=readmission,Frailty=TRUE,
	     hazard="Splines")

print(modSha,type="surv",conf=FALSE)

# /*** Cox proportional hazard model using Penalized likelihood ***/

modCox<-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
             n.knots=10,kappa1=10000,data=readmission,Frailty=FALSE,
	     hazard="Splines")
print(modCox)


}

\keyword{methods}
