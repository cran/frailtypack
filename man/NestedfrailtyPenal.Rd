

\name{frailtyPenal for Nested frailty models}
\alias{frailtyPenal for Nested frailty models}

\title{Fit Nested Frailty model using penalized likelihood estimation }

\description{Fit a nested frailty model using a Penalized Likelihood on the hazard function. Nested frailty models allow survival studies for hierarchically clustered data by including two iid gamma random effects.
Left truncated and censored data are allowed. Stratification analysis is allowed(maximum of strata=2).

The hazard function conditional on the two frailties \eqn{v_i} and \eqn{w_{ij}} for the \eqn{k^{th}} individual of the \eqn{j^{th}} subgroup of the \eqn{i^{th}} group is :

\deqn{\lambda_{ijk}(t|v_i,w_{ij})=v_iw_{ij}\lambda_0(t)exp(\bold{\beta^{'}X_{ijk}})}

\deqn{\small{ v_i\sim\Gamma\left(\frac{1}{\alpha},\frac{1}{\alpha}\right) \hspace{0.05cm}i.i.d. \hspace{0.2cm} \bold{E}(v_i)=1 \hspace{0.2cm}\bold{Var}(v_i)=\alpha \hspace{0.5cm}  w_{ij}\sim\Gamma\left(\frac{1}{\eta},\frac{1}{\eta}\right)\hspace{0.05cm}i.i.d. \hspace{0.2cm} \bold{E}(w_{ij})=1 \hspace{0.2cm}\bold{Var}(w_{ij})=\eta}}

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{X_{ijk}}
 denotes the covariate vector and \eqn{\beta} the corresponding vector of regression parameters.

}


\synopsis{

frailtyPenal(formula, formula.terminalEvent, data, Frailty = FALSE,
                 joint = FALSE, recurrentAG = FALSE, cross.validation =
                 FALSE, n.knots, kappa1, kappa2, maxit = 350)
}

\arguments{
  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right.  The response must be a
          survival object as returned by the 'Surv' function like in survival package. The \code{subcluster()} function is required.}

  \item{formula.terminalEvent}{Not required.}

  \item{data}{ a 'data.frame' in which to interpret the variables named in the
 'formula'.}

  \item{Frailty}{Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. If not, Cox proportional hazards model is estimated using Penalized Likelihood on the hazardfunction. The default is FALSE.}
 
  \item{joint}{Not required}

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
          second derivative of hazard function in the fit (penalized log likelihood). We advise the user 
          to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

  \item{kappa2}{ positive smoothing parameter for the second stratum, when data are stratified. See kappa1. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

}


\details{
  The estimated parameter are obtained by maximizing the penalized loglikelihood using the robust Marquardt algorithm (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm and a steepest descent algorithm.
 When frailty parameter is small, numerical problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see 
  Rondeau, 2003 for further details). 
  Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.
The smoothing parameter can be fixed a priori or chosen by maximizing a likelihood cross validation criterion.
 The iterations are stopped when the difference between two
 consecutive loglikelhoods was small \eqn{(<10^{-4})}, the estimated coefficients
 were stable (consecutive values \eqn{(<10^{-4})}, and the gradient small
 enough \eqn{(<10^{-6})}.
To be sure of having a positive function at all stages of
 the algorithm, the spline coefficients were reparametrizes to be positive
 at each stage. The variance space of the two random effects is reduced,
 so the variances are positive, and the correlation coefficient values are
 constrained to be  between -1 and 1.
The integrations in the full log likelihood werre evaluated using Gaussian quadrature. Laguerre polynomials with 10 points were used to treat the integrations on \eqn{[0,\infty[}.

\bold{INITIAL VALUES}

When a nested frailty model is fitted, at first, a Cox proportional hazards model with covariates (without random effects) is fitted and it provides initial values for the regression coefficients. Then, a model with covariates and only subgroup frailties is fitted, so it gives new initial values for the regression coefficients and for the subgroup frailty terms variance.  Then, a model with covariates and only group frailty terms is fitted, so it gives new initial values for the regression coefficients and for the group frailties variance. The last step is the maximisation of the penalized log-likelihood for all parameters.



\bold{PARAMETERS LIMIT VALUES}

As frailtypack is written in Fortran 77 some parameters had to be hard coded in. The default values of these parameters are, with the corresponding variable name in the fortran code between brackets.

maximum number of observations (ndatemax): 30000\cr
maximum number of groups (ngmax): 1000\cr
maximum number of subjects (nsujetmax): 15000\cr
maximum number of parameters (npmax) :50\cr
maximum number of covariates (nvarmax):50\cr
maximum number of subgroups (nssgmax):5000\cr


If these parameters are not large enough (an error message will let you know this), you need to reset them in nested.f and recompile. 
}

\value{
  a Nested frailty model or more generally an object of class 'frailtyPenal'. Methods defined for 'frailtyPenal' objects are provided for print, plot and summary. The following components are included in a 'frailtyPenal' object for nested frailty models.

\item{alpha}{variance of the cluster effect \eqn{(\bold{Var}(v_{i}))}}
\item{call}{The code used for fitting the model}
\item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
\item{cross.Val}{Logical value. Is cross validation procedure used for estimating the smoothing parameters?}
\item{DoF}{degrees of freedom}
\item{eta}{variance of the subcluster effect \eqn{(\bold{Var}(w_{ij}))}}
\item{formula}{the formula part of the code used for the model}
\item{groups}{the maximum number of groups used in the fit}
\item{subgroups}{the maximum number of subgroups used in the fit}
\item{kappa}{ A vector with the smoothing parameters corresponding to each baseline function as components}
\item{lam}{matrix of hazard estimates at x1 times and confidence bands. }
\item{lam2}{the same value as lam for  the second stratum}
\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit}
\item{n.iter}{number of iterations needed to converge}
\item{n.knots}{number of knots}
\item{n.strat}{A vector with the number of covariates of each type of hazard function as components}
\item{nvar}{A vector with the number of covariates of each type of hazard function as components}
\item{surv}{matrix of baseline survival estimates at x1 times and confidence bands. }
\item{surv2}{the same value as surv for the second stratum}
\item{type}{a character string specifying the type of censoring. Possible values are "right", "left", "counting", "interval", "interval2". The default is "right" or "counting" depending on whether the 'time2' argument is absent (not interval-censored data) or present (interval-censored data), respectively.}
\item{varH}{the variance matrix of all parameters (alpha, eta, the regression coefficients and the spline coefficients).}
\item{varHIH}{the robust estimation of the variance matrix of all parameters (alpha, eta, the regression coefficients and the spline coefficients).}
\item{x1}{vector of times where both survival and hazard function are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{the same value as x1 for the second stratum}







 }

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}


\references{

V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum penalized likelihood estimation.
\emph{Statistics in Medecine}, \bold{25}, 4036-4052.

V. Rondeau, D Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
\emph{Lifetime Data Analysis} \bold{9}, 139-153.


D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 


 }

 

\seealso{
 \code{\link{print.nestedPenal}},
 \code{\link{summary.nestedPenal}},
 \code{\link{plot.nestedPenal}},
 \code{\link{cluster}},
 \code{\link{subcluster}},
 \code{\link{strata}}

}


\examples{



 ### Nested model (or hierarchical model) with 2 covariates ###


\dontrun{

  data(dataNested)
  modClu<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+subcluster(subgroup)+
                  cov1+cov2,Frailty=TRUE,data=dataNested,n.knots=8,kappa1=50000)


# It takes around 24 minutes to converge (depends on the processor)#

  print(modClu)
  summary(modClu)
  plot(modClu)
  
  modClu.str<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+subcluster(subgroup)+
		  cov1+strata(cov2)
                  ,Frailty=TRUE,data=dataNested,n.knots=8,kappa1=20000,kappa2=20000)
  
 # It takes around 8 minutes to converge (depends on the processor)#
 
  print(modClu.str)
  summary(modClu.str)
  plot(modClu.str)
  
  
}


}

\keyword{models}

