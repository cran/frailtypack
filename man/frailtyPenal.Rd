
\name{frailtyPenal for Shared frailty models}
\alias{frailtyPenal}
\alias{frailtyPenal for Shared frailty models}

\title{ Fit a Shared Gamma Frailty model using a semi-parametrical penalized likelihood estimation or parametrical estimation }

\description{Fit a shared gamma frailty model using a semi-parametrical Penalized Likelihood estimation or parametrical estimation on the hazard function.
Left truncated and censored data and strata (max=2) are allowed.
It allows to obtain a non-parametrical smooth hazard of survival function.
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
            FALSE, n.knots, kappa1, kappa2, maxit = 350,hazard,
	    nb.int1,nb.int2)
}

\arguments{
  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right.  The response must be a
          survival object as returned by the 'Surv' function like in survival package.}

  \item{formula.terminalEvent}{Not required.}

  \item{data}{ a 'data.frame' in which to interpret the variables named in the
 'formula'.}
  
  \item{Frailty}{Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. If not, Cox proportional hazard model is estimated using Penalized Likelihood on the hazard function. The default is TRUE}
 
  \item{joint}{Not required.} 

  \item{recurrentAG}{Logical value. Is Andersen-Gill model fitted? If so indicates that recurrent event times with the 
                 counting process approach of Andersen and Gill is used. This formulation can be used for dealing with
                 time-dependent covariates. The default is FALSE.}

  \item{cross.validation}{Logical value. Is cross validation procedure used for estimating smoothing parameter 
                 in the penalized likelihood estimation?
                 If so a search of the smoothing parameter using cross validation is done, with kappa1 as the seed.
                 The cross validation is not implemented for two strata. The cross validation has been implemented 
                 for a Cox proportional hazard model, with no covariates. The default is FALSE.}

  \item{n.knots}{ integer giving the number of knots to use. Value required in the penalized likelihood estimation.
  It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions. 
  Number of knots must be between 4 and 20. (See Note)}

  \item{kappa1}{ positive smoothing parameter in the penalized likelihood estimation. The coefficient kappa of the integral of the squared
          second derivative of hazard function in the fit (penalized log likelihood). To obtain an initial 
value for \code{kappa1} (or \code{kappa2}), a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation). 
We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

  \item{kappa2}{positive smoothing parameter in the penalized likelihood estimation for the second stratum, when data are stratified. See kappa1. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}
  
 \item{hazard}{Type of hazard functions: "Splines" for semi-parametrical hazard functions with the penalized likelihood estimation, "Piecewise-per" for piecewise constant hazard function using percentile, "Piecewise-equi" for piecewise constant hazard function using equidistant intervals, "Weibull" for parametrical weibull functions. Default is "Splines".}
 
\item{nb.int1}{Number of intervals (between 1 and 20) for the parametrical hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nb.int2}{Argument used only for the Joint frailty model.}
}

\details{
  The estimated parameter are obtained using the robust Marquardt algorithm (Marquardt, 1963) which is a combination
  between a Newton-Raphson algorithm and a steepest descent algorithm. When frailty parameter is small, numerical
  problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see 
  Rondeau, 2003 for further details). 
  Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.1. The program fits, firstly, an adjusted Cox model to give new initial values for the splines and the regression coefficients. The variance of the frailty term \eqn{\theta} is initialized to 0.1. Then, a shared frailty model is fitted.

}

\value{

The following components are included in a 'frailtyPenal' object for shared frailty models.
\item{b}{sequence of the corresponding estimation of the coefficients for the hazard functions (parametrical or semi-parametrical), the random effects variances and the regression coefficients.}
\item{call}{The code used for the model.}
\item{coef}{the regression coefficients.}
\item{cross.Val }{Logical value. Is cross validation procedure used for estimating the smoothing parameters in the penalized likelihood estimation?  }
\item{DoF}{Degrees of freedom associated with the "kappa".}
\item{formula }{the formula part of the code used for the model.}
\item{groups}{the maximum number of groups used in the fit.}
\item{kappa}{ A vector with the smoothing parameters  in the penalized likelihood estimation corresponding to each baseline function as components.}
\item{lam}{matrix of hazard estimates and confidence bands.}
\item{lam2}{the same value as lam for the second stratum.}
\item{loglikPenal}{the complete marginal penalized log-likelihood in the semi-parametrical case.}
\item{loglik}{the marginal log-likelihood in the parametrical case.}
\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit.}
\item{n.iter}{number of iterations needed to converge.}
\item{n.knots}{number of knots for estimating the baseline functions in the penalized likelihood estimation.}
\item{n.strat}{ number of stratum.}
\item{surv}{matrix of baseline survival estimates and confidence bands.}
\item{surv2}{the same value as surv for the second stratum.}
\item{theta}{variance of frailty parameter.}
 \item{varH}{the variance matrix of all parameters before positivity constraint transformation (theta, the regression coefficients and the spline coefficients). Thenafter, the delta method is needed to obtain the estimated variance parameters.}
 \item{varHIH}{the robust estimation of the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
\item{x1}{vector of times where both survival and hazard function are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{vector of times for the second stratum (see x1 value).}
\item{type.of.hazard}{Type of hazard functions (0:"Splines", "1:Piecewise", "2:Weibull").}
\item{type.of.Piecewise}{Type of Piecewise hazard functions (1:"percentile", 0:"equidistant").}
\item{nbintervR}{Number of intervals (between 1 and 20) for the parametrical hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{npar}{number of parameters.}
\item{nvar}{number of explanatory variables.}
\item{noVar}{indicator of explanatory variables.}
\item{LCV}{the approximated likelihood cross-validation criterion in the semi-parametrical case (with H minus the converged hessien matrix, and l(.) the full log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{the Akaike information Criterion for the parametrical case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
\item{n.knots.temp}{initial value for the number of knots.}
\item{shape.weib}{shape parameter for the weibull hazard function.}
\item{scale.weib}{scale parameter for the weibull hazard function.}
\item{martingale.res}{martingale residuals for each cluster.}
\item{martingaleCox}{martingale residuals for observation in the Cox model.}
\item{Frailty}{Logical value. Was model with frailties fitted.}
\item{frailty.pred}{empirical Bayes prediction of the frailty term (ie, using conditional posterior distributions).}
\item{frailty.var}{variance of the empirical Bayes prediction of the frailty term.}
\item{frailty.sd}{standard error of the frailty empirical Bayes prediction.}
\item{linear.pred}{linear predictor: uses simply "Beta'X" in the cox proportional hazard model or "Beta'X + log z_i" in the shared Frailty models.}

\item{global_chisq}{a vector with the values of each multivariate Wald test.}
\item{dof_chisq}{a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test}{a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq}{a vector with the p_values for each global multivariate Wald test.}
\item{names.factor}{Names of the "as.factor" variables.}
}

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}



\references{
V. Rondeau, Y. Mazroui and J. R. Gonzalez. Frailtypack: An R package for the analysis of correlated survival data with frailty models using penalized likelihood estimation.
\emph{Journal of Statistical Software, In press}, 2012.

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
             n.knots=10,kappa1=10000,kappa2=10000,data=readmission,
	     Frailty=TRUE)

###  Shared model - recurrentAG=TRUE ###

frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
           as.factor(charlson)+cluster(id),data=readmission,Frailty=TRUE,
           n.knots=6,kappa1=100000,recurrentAG=TRUE)

###  Shared model - cross.validation=TRUE ###

\dontrun{ 
frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
          as.factor(charlson)+cluster(id),data=readmission, Frailty=TRUE,
          n.knots=6,kappa1=5000,recurrentAG=TRUE,cross.validation=TRUE)
}

}

\keyword{models}

