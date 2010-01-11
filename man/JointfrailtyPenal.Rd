
\name{frailtyPenal for Joint frailty models}
\alias{frailtyPenal for Joint frailty models}
\alias{frailtyPenal for Joint models}


\title{ Fit Joint Frailty model for recurrent and terminal events using
 penalized likelihood estimation }

\description{
Fit a joint frailty model for recurrent and terminal events using a Penalized Likelihood on the hazard function. Left-truncated and right-censored data are allowed. Stratified analysis is not possible. Joint frailty models allow studying, jointly, survival processes of recurrent and terminal events, by considering the terminal event 
as an informative censoring. A common frailty term \eqn{(\omega_i)}  for the two rates will take into account the heterogeneity in the data, associated with unobserved covariates. The frailty term acts differently for the two rates ( \eqn{\omega_i} for the recurrent rate and \eqn{\omega_i^{\alpha}} for death rate). The covariates could be different for the recurrent rate and death rate.


The joint model for the hazard functions for recurrent event \eqn{r_i(.)} and death \eqn{\lambda_i(.)} is :


\deqn{ r_i(t|\omega_i)=\omega_ir_0(t)exp(\bold{\beta_1^{'}Z_i(t)})}
 
\deqn{\lambda_i(t|\omega_i)=\omega_i^{\alpha}\lambda_0(t)exp(\bold{\beta_2^{'}Z_i(t)})}

where \eqn{r_0(t)} (resp. \eqn{\lambda_0(t)}) is the recurrent (resp. terminal) event baseline hazard function, \eqn{\bold{\beta_1}} (resp. \eqn{\bold{\beta_2}}) the regression coefficient vector, \eqn{\bold{Z_i(t)}} the covariate vector. And where \eqn{\omega_i\sim\bold{\Gamma}(\frac{1}{\theta},\frac{1}{\theta})}  is iid. 

}


\synopsis{

frailtyPenal(formula, formula.terminalEvent, data, Frailty = FALSE,
                 joint = FALSE, recurrentAG = FALSE, cross.validation =
                 FALSE, n.knots, kappa1, kappa2, maxit = 350)
}

\arguments{
  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right.  The response must be a
          survival object as returned by the 'Surv' function like in survival package.}



  \item{formula.terminalEvent}{ a formula object, only requires terms on the right to indicate which variables are modelling the terminal event.  }

  \item{data}{ a 'data.frame' in which to interpret the variables named in the
 'formula' and 'formula.terminalEvent'.}


  \item{Frailty}{Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. 
                 If not, Cox proportional hazards model is estimated using Penalized Likelihood on the hazard
                 function. The default is FALSE}
 
  \item{joint}{Logical value. Is joint model fitted? If so, 'formula.terminalEvent' is required. The default is FALSE}

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
          second derivative of hazard function in the fit (penalized log likelihood). To obtain a good
value for \code{kappa1}, a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation) with the event indicator as the event of interest. 
We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

  \item{kappa2}{ positive smoothing parameter for the terminal event rate. To obtain a good
value for \code{kappa2}, a solution is to fit the corresponding Cox model using cross-validation (See cross.validation) with the death indicator as the event of interest. See kappa1. Stratification is not allowed here. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

}




\details{

  The estimated parameter are obtained by maximizing the penalized loglikelihood  using the robust Marquardt algorithm (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm and a steepest descent algorithm. When frailty parameter is small, numerical problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see Rondeau, 2003 for further details). 
  Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.
The inverse of the hessian matrix is the variance estimator and to deal with thepositivity constraint of the variance component, a squared transformation is used and the Standard errors are computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth parameter can be chosen by maximizing a likelihood cross validation criterion (Joly and other, 1998). The integrations in the full log likelihood werre evaluated using Gaussian quadrature. Laguerre polynomials with 20 points were used to treat the integrations on \eqn{[0,\infty[}.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.5. The program fits an adjusted Cox model to have new initial values for the regression and the splines coefficients. The variance of the frailty term \eqn{\theta} and the coefficient \eqn{\alpha} associated in the death hazard function are initialized to 1. Then, it fits a joint frailty model.


\bold{PARAMETERS LIMIT VALUES}

As frailtypack is written in Fortran 77 some parameters had to be hard coded in.The default values of these parameters are, with the corresponding variable namein the fortran code between brackets.


maximum number of observations(ndatemax): 30000\cr
maximum number of groups(ngmax): 1000\cr
maximum number of subjects(nsujetmax): 15000\cr
maximum number of parameters (npmax) :50\cr
maximum number of covariates (nvarmax):50\cr

If these parameters are not large enough (an error message will
 let you know this), you need to reset them in joint.f and recompile. 
}

\value{
Parameters estimates of a joint frailty model, more generally a 'fraityPenal' object. Methods defined for 'frailtyPenal' objects are
provided for print,  plot and summary. The following components are included in a 'frailtyPenal' object for Joint frailty models.

\item{alpha}{the coefficient associated with the frailty parameter in the terminal hazard function}
\item{call}{The code used for fitting the model}
\item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }

\item{cross.Val}{Logical value. Is cross validation procedure used for estimating the smoothing parameters?}

\item{formula}{the formula part of the code used for the model}
\item{groups}{the maximum number of groups used in the fit}
\item{kappa}{ A vector with the smoothing parameters corresponding to each baseline function as components}
\item{lam}{matrix of hazard estimates at x1 times and confidence bands. }
\item{lam2}{the same value as lam for  the second stratum}

\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit}
\item{n.iter}{number of iterations needed to converge}
\item{nvar}{A vector with the number of covariates of each type of hazard function as components}
\item{surv}{matrix of baseline survival estimates at x1 times and confidence bands. }
\item{surv2}{the same value as surv for the  the second stratum}
\item{theta}{variance of the frailty parameter \eqn{(\bold{Var}(\omega_i))}}
\item{type}{character string specifying the type of censoring. Possible values are "right", "left", "counting", "interval", "interval2". The default is "right" or "counting" depending on whether the 'time2' argument is absent (not interval-censored data) or present (interval-censored data), respectively.}
 \item{varH}{the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
 \item{varHIH}{the robust estimation of the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
\item{x1}{vector of times where both survival and hazard function for the recurrent events are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{see x1 value for the terminal events survival and hazard function}

}

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}


\references{

V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P. Soubeyran (2007). Joint frailty models for recurring events and death using maximum penalized likelihood estimation:application on cancer events.
\emph{Biostatistics}, \bold{8},4, 708-721.

V. Rondeau, D Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
\emph{Lifetime Data Analysis} \bold{9}, 139-153.

D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 

 }
 
\seealso{
 \code{\link{summary.jointPenal}},
 \code{\link{print.jointPenal}},
 \code{\link{plot.jointPenal}}, 
 \code{\link{readmission}},
 \code{\link{terminal}},
 \code{\link{cluster}}

 }


\examples{

### Joint model (recurrent and terminal events) with 2 covariates ###
### on a simulated dataset ###


\dontrun{

  data(readmission)

## Gap-time ##
modJoint_gap<-frailtyPenal(Surv(time,event)~cluster(id)+sex+as.factor(dukes)
		+as.factor(charlson)+terminal(death),
		formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
		data=readmission,n.knots=14,kappa1=9550000000,
		kappa2=1410000000000,Frailty=TRUE,joint=TRUE,recurrentAG=FALSE)

## Calendar time ##
modJoint_calendar<-frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
		+as.factor(dukes)+as.factor(charlson)+terminal(death),
		formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
		data=readmission,n.knots=10,kappa1=9550000000,
		kappa2=1410000000000,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE)
    
  print(modJoint_gap)
  summary(modJoint_gap)
  plot(modJoint_gap)
  
  print(modJoint_calendar)
  summary(modJoint_calendar)
  plot(modJoint_calendar)

}

# A model takes around 1 minute to converge #


}

\keyword{models}
