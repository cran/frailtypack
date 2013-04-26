\name{print.multivePenal}
\Rdversion{1.1}
\alias{print.multivePenal}

\title{
Print a Short Summary of parameter estimates of a multivariate frailty model
}
\description{
  Prints a short summary of parameter estimates of a multivariate frailty model, or more generally an object of class 'multivePenal'.
}
\usage{
\method{print}{multivePenal}(x, digits = max(options()$digits - 4, 6), ...)
}

\arguments{
  \item{x}{the result of a call to the multivePenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}

\value{
Print, separately for each type of event (recurrent1, recurrent2 and terminal), the parameter estimates of the survival or hazard functions.

 \item{n}{the number of observations used in the fit.}
 \item{n.groups}{the maximum number of groups used in the fit}
 \item{n.events}{the number of recurrent event of type 1 observed in the fit}
 \item{n.deaths}{the number of terminal events observed in the fit}
 \item{n.events2}{the number of the recurrent event of type 2 observed in the fit}
 \item{theta}{variance of the common frailty term}
 \item{alpha}{the coefficient associated with the frailty parameter terminal event hazard function}
 \item{coef}{the coefficients of the linear predictor, which multiply the columns of the model matrix. }
 \item{varH}{the variance matrix of the dependency parameters (theta, eta, alpha1, alpha2) and of the coefficients.}
 \item{varHIH}{the robust estimation of the dependency parameters and of the coefficients.}
 \item{SE(H)}{the standard error of the estimates deduced from the variance matrix of the dependency parameters and of the coefficients.}
 \item{SE(HIH)}{the standard error of the estimates deduced from the robust estimation of the variance matrix of the dependency parameters and of the coefficients.}
 \item{p}{p-value}
}

\seealso{
 \code{\link{multivePenal}}
}


\examples{
# /*** multivariate frailty model ***/
\dontrun{

data(dataMultiv)
dataMultiv$timegap<-dataMultiv$t1-dataMultiv$t0

## Calendar-time ##


MultivSpli <-multivePenal(Surv(timegap,deltaEvent1)~cluster(id)+var1+var2+terminal(deltadc)+event2(deltaEvent2),formula.terminalEvent=~var1,formula2=~var1+var2+var3,data=dataMultiv,Frailty=TRUE, recurrentAG=F,cross.validation=T,n.knots=c(4,4,4),kappa=c(1,1,1),maxit=20,hazard="Splines",initialize=T)

## print a fit
MultivSpli

}
}

\keyword{methods}
\keyword{multive}
