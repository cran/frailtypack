
\name{print.additivePenal}
\Rdversion{1.1}

\alias{print.additivePenal}

\title{Print a Short Summary of parameter estimates of an additive frailty model}
\description{
  Prints a short summary of the parameter estimates of an additive frailty model or more generally of an 'additivePenal' object 
}
\usage{
\method{print}{additivePenal}(x, digits = max(options()$digits - 4, 6), ...)
}
\arguments{
  \item{x}{the result of a call to the additivePenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\value{

 \item{n}{the number of observations used in the fit.}
 \item{n.groups}{the maximum number of groups used in the fit}
 \item{n.events}{the number of events observed in the fit}
 \item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
 \item{SE(H)}{the standard error deduced from the variance matrix of theta and of the coefficients.}
 \item{SE(HIH)}{the standard error deduced from the robust estimation of the variance matrix of theta and of the coefficients.}
 \item{z}{quantile of the \eqn{\chi^2} mixture distribution }
 \item{p}{p-value}
 \item{Variance for the random intercept}{Variance for the random effect associated to the baseline risk functon}
 \item{Variance for the random slope}{Variance for the random effect associated to the treatment effects across trials}


}


\seealso{
 \code{\link{additivePenal}}
}

\examples{

# Additive frailty model with 1 covariate

\dontrun{

data(dataAdditive)
modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=862,hazard="Splines")

# 'var1' is boolean as a treatment variable.

print(modAdd)

}

}

\keyword{methods}
