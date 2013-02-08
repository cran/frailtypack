
\name{print.jointPenal}
\Rdversion{1.1}

\alias{print.jointPenal}

\title{Print a Short Summary of parameter estimates of a joint frailty model}
\description{
  Prints a short summary of parameter estimates of a joint frailty model, or more generally an object of class 'frailtyPenal' for joint frailty models. 
}
\usage{

\method{print}{jointPenal}(x, digits = max(options()$digits - 4, 6), ...)

}
\arguments{
  \item{x}{the result of a call to the jointPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\value{

Print, separately for each type of event (recurrent and terminal), the parameter estimates of the survival or hazard functions.

 \item{n}{the number of observations used in the fit.}
 \item{n.groups}{the maximum number of groups used in the fit}
 \item{n.events}{the number of events observed in the fit}

 \item{theta}{variance of the common frailty term}
 \item{alpha}{the coefficient associated with the frailty parameter terminal event hazard function}
 \item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
 \item{varH}{the variance matrix of theta and of the coefficients.}
 \item{varHIH}{the robust estimation of the variance matrix of theta and of the coefficients.}
 \item{SE(H)}{the standard error of the estimaes deduced from the variance matrix of theta and of the coefficients.}
 \item{SE(HIH)}{the standard error of the estimaes deduced from the robust estimation of the variance matrix of theta and of the coefficients.}
 \item{p}{p-value}


}


\seealso{ \code{\link{summary.jointPenal}},  \code{\link{frailtyPenal for Joint frailty models}},\code{\link{plot.jointPenal}} }

\examples{

# /*** Joint frailty model ***/

\dontrun{

data(readmission)

## Gap-time ##
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes
                +charlson+terminal(death),formula.terminalEvent=~sex
                +dukes+charlson,data=readmission,n.knots=14,
                kappa1=9.55e+9,kappa2=1.41e+12,Frailty=TRUE,joint=TRUE)

# It takes around 1 minute to converge #

print(modJoint.gap)

}

}

\keyword{methods}
