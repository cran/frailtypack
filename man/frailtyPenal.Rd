\name{frailtyPenal}
\alias{frailtyPenal}
\title{ Fit Shared Gamma Frailty model using penalized likelihood estimation }
\description{Fit a shared gamma frailty model using a Penalized Likelihood on the hazard function.
Left truncated and censored data and strata (max=2) are allowed.
It allows to obtain a non-parametric smooth hazard of survival function.
This approach is different from the partial penalized likelihood approach of Therneau et al. 
 
}
\usage{
frailtyPenal(formula, data, Frailty = TRUE, n.knots, 
               kappa1, kappa2, maxit=350)
}

\arguments{
  \item{formula}{ a formula object, with the response on the left of a '~'
          operator, and the terms on the right.  The response must be a
          survival object as returned by the 'Surv' function like in survival package. 
}
  \item{data}{ a data.frame in which to interpret the variables named in the 'formula'.}
  \item{Frailty}{Logical value. Is model with frailties fitted? If so variance of frailty parameter is estimated. 
                 If not, Cox proportional hazards model is estimated using Penalized Likelihood on the hazard
                 function}
  \item{n.knots}{ integer giving the number of knots to use. Value required.
  It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions. 
  Number of knots must be between 4 and 20.}
  \item{kappa1}{ positive smoothing parameter. The coefficient kappa of the integral of the squared
          second derivative of hazard function in the fit (penalized log likelihood). We advise the user 
          to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.}
  \item{kappa2}{ positive smoothing parameter for the second stratum, when data are stratified. See kappa1. }
  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

}
\details{
  The estimated parameter are obtained using the robust Marquardt algorithm (Marquardt, 1963) which is a combination
  between a Newton-Raphson algorithm and a steepest descent algorithm. When frailty parameter is small, numerical
  problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see 
  Rondeau, 2003 for further details). 
Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.
}
\value{
  an object of class '"frailtyPenal"'. Methods defined for 'frailtyPenal' objects are
     provided for print and plot. The following components are included in a 'frailtyPenal'object.

  \item{n}{the number of observations used in the fit.}
  \item{groups}{the maximum number of groups used in the fit}
  \item{n.events}{the number of events observed in the fit}
  \item{logVerComPenal}{the complete marginal penalized log-likelihood }
  \item{theta}{variance of frailty parameter}
  \item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
  \item{varH}{the variance matrix of theta and of the coefficients.}
  \item{varHIH}{the robust estimation of the variance matrix of theta and of the coefficients.}
  \item{x1}{vector of times where both survival and hazard function are estimated. 
           By default seq(0,max(time),length=99), where time is the vector of survival times.}
  \item{lam}{matrix of hazard estimates at x1 times and confidence bands. }
  \item{surv}{matrix of baseline survival estimates at x1 times and confidence bands. }
  \item{x2}{see x1 value for the second stratum}
  \item{lam2}{the same value as lam for the second stratum}
  \item{surv2}{the same value as surv for the second stratum}
}



\references{D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 


V. Rondeau, D Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
\emph{Lifetime Data Analysis} \bold{9}, 139-153.

McGilchrist CA and Aisbett CW (1991). Regression with frailty in survival analysis.
\emph{Biometrics} \bold{47}, 461-466.

 }

 

\seealso{ \code{\link{print.frailtyPenal}}, \code{\link{summary.frailtyPenal}}}


\examples{
data(kidney)
#Shared frailty model
frailtyPenal(Surv(time,status)~sex+age+cluster(id),
             n.knots=12,kappa1=1000,data=kidney)


#model without frailties (e.g., Cox proportional hazards 
#                          estimated via penalized likelihood)
frailtyPenal(Surv(time,status)~sex+age+cluster(id),
             n.knots=12,kappa1=1000,data=kidney,Frailty=FALSE)


# truncated data

# first, we create a hypothetical truncated data
kidney$tt0<-rep(0,nrow(kidney))
kidney$tt0[1:3]<-c(2,9,13)

# then, we fit the model
frailtyPenal(Surv(tt0,time,status)~sex+age+cluster(id),
             n.knots=12,kappa1=1000,data=kidney)


#stratified data. Let's use another dataset
data(readmission)
frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
             n.knots=10,kappa1=10000,kappa2=10000,data=readmission)



}

\keyword{ survival }

