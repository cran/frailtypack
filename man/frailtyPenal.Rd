
\name{frailtyPenal for Shared frailty models}
\alias{frailtyPenal}
\alias{frailtyPenal for Shared frailty models}

\title{ Fit Shared Gamma Frailty model using penalized likelihood estimation }

\description{Fit a shared gamma frailty model using a Penalized Likelihood on the hazard function.
Left truncated and censored data and strata (max=2) are allowed.
It allows to obtain a non-parametric smooth hazard of survival function.
This approach is different from the partial penalized likelihood approach of Therneau et al. 
 
The hazard function, conditionnal on the frailty term \eqn{z_i}, of a shared frailty model for the \eqn{j^{th}} subject in the \eqn{i^{th}} group:

\deqn{\lambda_{ij}(t|z_i)=\lambda_0(t)z_iexp(\bold{\beta^{'}X_{ij}})}

\deqn{z_i\sim\Gamma\left(\frac{1}{\theta},\frac{1}{\theta}\right)
 \hspace{0.5cm} \bold{E}(z_i)=1 \hspace{0.5cm}\bold{Var}(z_i)=\theta}

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{\bold{\beta}} the vector of the regression coefficient associated to the covariate vector \eqn{\bold{X_{ij}}} for the \eqn{j^{th}} individual in the \eqn{i^{th}} group.



}

\usage{

frailtyPenal(formula, formula.terminalEvent, data, Frailty = FALSE,
                 joint = FALSE, recurrentAG = FALSE, cross.validation =
                 FALSE, n.knots, kappa1, kappa2, maxit = 350)
}

\arguments{
  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right.  The response must be a
          survival object as returned by the 'Surv' function like in survival package.}

  \item{formula.terminalEvent}{Not required.}

  \item{data}{ a 'data.frame' in which to interpret the variables named in the
 'formula'.}
  
  \item{Frailty}{Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. If not, Cox proportional hazards model is estimated using Penalized Likelihood on the hazard function. The default is TRUE}
 
  \item{joint}{Not required.} 

  \item{recurrentAG}{Logical value. Is Andersen-Gill model fitted? If so indicates that recurrent event times with the 
                 counting process approach of Andersen and Gill is used. This formulation can be used for dealing with
                 time-dependent covariates. The default is FALSE.}

  \item{cross.validation}{Logical value. Is cross validation procedure used for estimating smoothing parameter?
                 If so a search of the smoothing parameter using cross validation is done, with kappa1 as the seed.
                 The cross validation is not implemented for two strata. The cross validation has been implemented 
                 for a Cox proportional hazard model, with no covariates. The default is FALSE.}

  \item{n.knots}{ integer giving the number of knots to use. Value required.
  It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions. 
  Number of knots must be between 4 and 20.(See Note)}

  \item{kappa1}{ positive smoothing parameter. The coefficient kappa of the integral of the squared
          second derivative of hazard function in the fit (penalized log likelihood). To obtain an initial 
value for \code{kappa1} (or \code{kappa2}), a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation). 
We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

  \item{kappa2}{positive smoothing parameter for the second stratum, when data are stratified. See kappa1. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

}

\details{
  The estimated parameter are obtained using the robust Marquardt algorithm (Marquardt, 1963) which is a combination
  between a Newton-Raphson algorithm and a steepest descent algorithm. When frailty parameter is small, numerical
  problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see 
  Rondeau, 2003 for further details). 
  Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.1. The program fits, firstly, an adjusted Cox model to give new initial values for the splines and the regression coefficients. The variance of the frailty term \eqn{\theta} is initialized to 0.1. Then, a shared frailty model is fitted.


\bold{PARAMETERS LIMIT VALUES}

As frailtypack is written in Fortran 77 some parameters had to be hard coded in. The default values of these parameters are, with the corresponding variable name in the fortran code between brackets.

maximum number of observations (ndatemax): 30000\cr
maximum number of groups (ngmax): 5000\cr
maximum number of subjects (nsujetmax): 20000\cr
maximum number of parameters (npmax) :50\cr
maximum number of covariates (nvarmax):50\cr



If these parameters are not large enough (an error message will let you know this), you need to reset
them in frailtypack.f and recompile. In particular, the statements defining these parameters are


}

\value{

The following components are included in a 'frailtyPenal' object for shared frailty models.

\item{call}{The code used for the model}
\item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
\item{cross.Val }{Logical value. Is cross validation procedure used for estimating the smoothing parameters?  }
\item{DoF }{Degrees of freedom associated with the "kappa"}
\item{formula }{the formula part of the code used for the model  }
\item{groups}{the maximum number of groups used in the fit}
\item{kappa  }{ A vector with the smoothing parameters corresponding to each baseline function as components }
\item{lam}{matrix of hazard estimates at x1 times and confidence bands. }
\item{lam2}{the same value as lam for the second stratum}
\item{loglikPenal}{the complete marginal penalized log-likelihood }
\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit}
\item{n.iter  }{number of iterations needed to converge  }
\item{n.knots }{number of knots for estimating the baseline function  }
\item{n.strat  }{ number of stratum }
\item{surv}{matrix of baseline survival estimates at x1 times and confidence bands. }
\item{surv2}{the same value as surv for the second stratum}
\item{theta}{variance of frailty parameter}
\item{type}{charater string specifying the type of censoring. Possible values are "right", "left", "counting", "interval", "interval2". The default is "right" or "counting" depending on whether the 'time2' argument is absent or present, respectively.}
 \item{varH}{the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
 \item{varHIH}{the robust estimation of the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
\item{x1}{vector of times where both survival and hazard function are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{see x1 value for the second stratum}

}

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}



\references{

V. Rondeau, D Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
\emph{Lifetime Data Analysis} \bold{9}, 139-153.

McGilchrist CA and Aisbett CW (1991). Regression with frailty in survival analysis.
\emph{Biometrics} \bold{47}, 461-466.

D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 

 }

 

\seealso{
\code{\link{print.frailtyPenal}},
\code{\link{summary.frailtyPenal}},
\code{\link{plot.frailtyPenal}},
\code{\link{cluster}},
\code{\link{strata}},
\code{\link[survival]{Surv}}
}


\examples{


 ###  Shared model  ###


data(kidney)
frailtyPenal(Surv(time,status)~ cluster(id)+sex+age,
             n.knots=12,kappa1=10000,data=kidney,Frailty=TRUE)

 ###  COX proportional hazard model (SHARED without frailties) ### 
 ### estimated with penalized likelihood ###

frailtyPenal(Surv(time,status)~sex+age+cluster(id),
             n.knots=12,kappa1=10000,data=kidney,Frailty=FALSE)

###  Shared model with truncated data ###

# Here is created a hypothetical truncated data
kidney$tt0<-rep(0,nrow(kidney))
kidney$tt0[1:3]<-c(2,9,13)

# then, we fit the model
frailtyPenal(Surv(tt0,time,status)~sex+age+cluster(id),
             n.knots=12,kappa1=10000,data=kidney,Frailty=TRUE)

###  Shared model - stratified analysis ###

data(readmission)
frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
             n.knots=10,kappa1=10000,kappa2=10000,data=readmission,Frailty=TRUE)

###  Shared model - recurrentAG=TRUE ###

frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
             as.factor(charlson)+cluster(id),data=readmission, Frailty=TRUE,
             n.knots=6,kappa1=100000,recurrentAG=TRUE)

###  Shared model - cross.validation=TRUE ###

\dontrun{ 
frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
             as.factor(charlson)+cluster(id),data=readmission, Frailty=TRUE,
             n.knots=6,kappa1=5000,recurrentAG=TRUE,cross.validation=TRUE)
}

}

\keyword{models}

