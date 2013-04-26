\name{multivePenal}
\alias{multivePenal}
\alias{transfo.table}
\alias{multivePenal for multivariate frailty model}

\title{
Fit a multivariate frailty model for two types of recurrent events and a terminal event using semiparametric penalized likelihood estimation or a parametrical estimation.
}
\description{
Fit a multivariate frailty model for two types of recurrent events with a terminal event using a penalized likelihood estimation on the hazard function or a parametric estimation. Right-censored data is allowed. Left-truncated data and stratified analysis are not possible.
Multivariate frailty models allow studying, jointly, three survival dependent processes for two types recurrent events and a terminal event, by considering the terminal event as an informative censoring. 
Dependencies between these three types of event are taken into account by \eqn{\alpha_1,\alpha_2,\eta,\theta}. If \eqn{\alpha_1} and \eqn{\theta} are both significantly different from 0, then the recurrent events of type 1 and death are significantly associated (the sign of the association is the sign of \eqn{\alpha_1}).
If \eqn{\alpha_2} and \eqn{\eta} are both significantly different from 0, then the recurrent events of type 2 and death are significantly associated (the sign of the association is the sign of \eqn{\alpha_2}).
If \eqn{rho} is significantly different from 0, then the recurrent events of type 1 and the recurrent events of type 2 are significantly associated (the sign of the association is the sign of \eqn{rho}).


The multivariate frailty model for two types of recurrent events with a terminal event is (in the calendar or time-to-event timescale):

\deqn{\left\{
\begin{array}{lll}
r_{i}^{(1)}(t|u_i,v_i) &= r_0^{(1)}(t)\exp({{\beta_1^{'}}}Z_{i}(t)+u_i) &\quad \mbox{(rec. of type 1)}\\
r_{i}^{(2)}(t|u_i,v_i) &= r_0^{(2)}(t)\exp({{\beta_2^{'}}}Z_{i}(t)+v_i) &\quad \mbox{(rec. of type 2)}\\
\lambda_i(t|u_i,v_i) &= \lambda_0(t)\exp({{\beta_3^{'}}}Z_{i}(t)+\alpha_1u_i+\alpha_2v_i) &\quad \mbox{(death)}\\
\end{array}
\right.
}

where \eqn{r_0^{(l)}(t)}, \eqn{l\in{1,2}} and \eqn{\lambda_0(t)} are respectively the recurrent and terminal event baseline hazard functions, and \eqn{\beta_1,\beta_2,\beta_3} the regression coefficient vectors associated with \eqn{Z_{i}(t)} the covariate vector. The covariates could be different for the different event hazard functions and may be time-dependent. We consider that death stops new occurrences of recurrent events of any type, hence given \eqn{t>D}, \eqn{dN^{R(l)*}(t), l\in{1,2}} takes the value 0. Thus, the terminal and the two recurrent event processes are not independent or even conditional upon frailties and covariates. We consider the hazard functions of recurrent events among individuals still alive. 
%
The three components in the above multivariate frailty model are linked together by two Gaussian and correlated random effects \eqn{u_i,v_i}:\\
%

\eqn{(u_i,v_i)^{T}\sim\mathcal{N}\left({{0}},\Sigma_{uv}\right)}, with \deqn{\Sigma_{uv}=\left(\begin{array}{cc}
\theta & \rho\sqrt{\theta\eta} \\
\rho\sqrt{\theta\eta}&\eta
\end{array}\right)}

}

\usage{

multivePenal(formula, formula.terminalEvent, formula2, data, Frailty = TRUE,
            initialize = TRUE, recurrentAG = FALSE, cross.validation = FALSE,
            n.knots, kappa, maxit = 350, hazard = "Splines", nb.int,
            print.times = T)

}

\arguments{
  \item{formula}{
 a formula object, with the response on the left of a \eqn{\texttildelow} operator, and the terms on the right. The response must be a survival object as returned by the 'Surv' function like in survival package.
}
  \item{formula.terminalEvent}{
a formula object, only requires terms on the right to indicate which variables are modelling the terminal event.
}
  \item{formula2}{
a formula object, only requires terms on the right to indicate which variables are modelling the second recurrent event.
}
  \item{data}{
 a 'data.frame' in which to interpret the variables named in the 'formula', 'formula.terminalEvent' and formula2.
}
  \item{Frailty}{
Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. The default is FALSE.
}
  \item{initialize}{
Logical value to initialize parameters. Default equals to 1.
}
  \item{recurrentAG}{
Logical value. Is Andersen-Gill model fitted? If so indicates that recurrent event times with the counting process approach of Andersen and Gill is used. This formulation can be used for dealing with time-dependent covariates. The default is FALSE.
}
  \item{cross.validation}{
Logical value. Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation? If so a search of the smoothing parameter using cross validation is done, with kappa1 as the seed. The cross validation is not implemented for two strata. The cross validation has been implemented for a Cox proportional hazard model, with no covariates. The default is FALSE.
}
  \item{n.knots}{
integer vector of length 3 giving the number of knots to use. Value required in the penalized likelihood estimation.
  It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions.
  Number of knots must be between 4 and 20. (See Note)
}
  \item{kappa}{
vector of length 3 for positive smoothing parameter in the penalized likelihood estimation. The coefficient kappa of the integral of the squared second derivative of hazard function in the fit (penalized log likelihood). To obtain a good value for \code{kappa1}, a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation) with the event indicator as the event of interest. We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity of the results to varying them. Value required.(See Note)
}
  \item{maxit}{
maximum number of iterations for the Marquardt algorithm. Default is 350
}
  \item{hazard}{
Type of hazard functions: "Splines" for semi-parametrical hazard functions with the penalized likelihood estimation, "Piecewise-per" for piecewise constant hazard function using percentile, "Piecewise-equi" for piecewise constant hazard function using equidistant intervals, "Weibull" for parametrical weibull function. Default is "Splines".
}
  \item{nb.int}{
An integer vector of length 3.first is the Number of intervals (between 1 and 20) for the recurrent parametrical hazard functions ("Piecewise-per", "Piecewise-equi"). Second is the Number of intervals (between 1 and 20) for the death parametrical hazard functions ("Piecewise-per", "Piecewise-equi"). Third is Number of intervals (between 1 and 20) for the second recurrent parametrical hazard functions ("Piecewise-per", "Piecewise-equi")
}
  \item{print.times}{
a logical parameter to print iteration process. Default is FALSE.
}
}

\value{
Parameters estimates of a multive joint frailty model, more generally a 'fraityPenal' object. Methods defined for 'frailtyPenal' objects are provided for print,  plot and summary. The following components are included in a 'multivePenal' object for multive Joint frailty models.

\item{b}{sequence of the corresponding estimation of the splines coefficients, the random effects variances, the power of the frailties and the regression coefficients.}
\item{call}{The code used for fitting the model.}
\item{n}{the number of observations used in the fit.}
\item{groups}{the maximum number of groups used in the fit.}
\item{n.events}{the number of recurrent events observed in the fit.}
\item{n.deaths}{the number of deaths events observed in the fit.}
\item{n.events2}{the number of the second recurrent events observed in the fit.}
\item{loglikPenal}{the complete marginal penalized log-likelihood in the semi-parametrical case.}
\item{loglik}{the marginal log-likelihood in the parametrical case.}
\item{LCV}{the approximated likelihood cross-validation criterion in the semi parametrical case (with H minus the converged hessien matrix, and l(.) the full log-likelihood.\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{the Akaike information Criterion for the parametrical case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
\item{theta}{variance of the frailty parameter \eqn{(\bold{Var}(u_i))}}
\item{eta}{variance of the frailty parameter \eqn{(\bold{Var}(v_i))}}
\item{alpha1}{the coefficient associated with the frailty parameter \eqn{u_i} in the terminal hazard function.}
\item{alpha2}{the coefficient associated with the frailty parameter \eqn{v_i} in the terminal hazard function.}
\item{rho}{the correlation coefficient between \eqn{u_i} and \eqn{v_i}}
\item{npar}{number of parameters.}
\item{coef}{the regression coefficients.}
\item{nvar}{A vector with the number of covariates of each type of hazard function as components.}
\item{varH}{the variance matrix of all parameters before positivity constraint transformation (theta, the regression coefficients and the spline coefficients). Thenafter, the delta method is needed to obtain the estimated variance parameters.}
\item{varHIH}{the robust estimation of the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
\item{formula}{the formula part of the code used for the model for the recurrent event.}
\item{formula.terminalEvent}{the formula part of the code used for the model for the terminal event.}
\item{formula2}{the formula part of the code used for the model for the second recurrent event.}
\item{x1}{vector of times for hazard functions of the recurrent events of type 1 are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{lam}{matrix of hazard estimates and confidence bands.}
\item{surv}{matrix of baseline survival estimates and confidence bands.}
\item{x2}{vector of times for the terminal event (see x1 value).}
\item{x3}{vector of times for the recurrent event of type 2 (see x1 value).}
\item{lam2}{the same value as lam for the terminal event.}
\item{surv2}{the same value as surv for the terminal event.}
\item{lam3}{the same value as lam for the recurrent event of type 2.}
\item{surv3}{the same value as surv for the recurrent event of type 2.}
\item{xSu1}{vector of times for the survival function of the recurrent event of type 1}
\item{xSu2}{vector of times for for the survival function of the terminal event}
\item{xSu3}{vector of times for the survival function of the recurrent event of type 2}
\item{type.of.Piecewise}{Type of Piecewise hazard functions (1:"percentile", 0:"equidistant").}
\item{n.iter}{number of iterations needed to converge.}
\item{type.of.hazard}{Type of hazard functions (0:"Splines", "1:Piecewise", "2:Weibull").}
\item{n.knots}{number of knots for estimating the baseline functions.}
\item{kappa}{ A vector with the smoothing parameters  in the penalized likelihood estimation corresponding to each baseline function as components.}
\item{cross.Val}{Logical value. Is cross validation procedure used for estimating the smoothing parameters  in the penalized likelihood estimation?}
\item{n.knots.temp}{initial value for the number of knots.}
\item{zi}{splines knots.}
\item{time}{knots for Piecewise hazard function for the recurrent event of type 1.}
\item{timedc}{knots for Piecewise hazard function for the terminal event.}
\item{time2}{knots for Piecewise hazard function for the recurrent event of type 2.}
\item{noVar}{indicator vector for reccurrent, death and recurrent 2 explanatory variables.}
\item{nvarRec}{number of the recurrent of type 1 explanatory variables.}
\item{nvarEnd}{number of death explanatory variables.}
\item{nvarRec2}{number of the recurrent of type 2 explanatory variables.}
\item{nbintervR}{Number of intervals (between 1 and 20) for the the recurrent of type 1 parametrical hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nbintervDC}{Number of intervals (between 1 and 20) for the death parametrical hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nbintervR2}{Number of intervals (between 1 and 20) for the the recurrent of type 2 parametrical hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{istop}{Vector of the convergence criteria.}
\item{shape.weib}{shape parameters for the weibull hazard function.}
\item{scale.weib}{scale parameters for the weibull hazard function.}

\item{martingale.res}{martingale residuals for each cluster (recurrent of type 1).}
\item{martingaledeath.res}{martingale residuals for each cluster (death).}
\item{martingale2.res}{martingale residuals for each cluster (recurrent of type 2).}
\item{frailty.pred}{empirical Bayes prediction of the first frailty term.}
\item{frailty2.pred}{empirical Bayes prediction of the second frailty term.}
\item{frailty.var}{variance of the empirical Bayes prediction of the first frailty term.}
\item{frailty2.var}{variance of the empirical Bayes prediction of the second frailty term.}
\item{frailty.corr}{Correlation between the empirical Bayes prediction of the two frailty.}
\item{linear.pred}{linear predictor: uses Beta'X + ui in the multivariate frailty models.}
\item{lineardeath.pred}{linear predictor for the terminal part form the multivariate frailty models: Beta'X + alpha1 ui + alpha2 vi}
\item{linear2.pred}{linear predictor: uses Beta'X + vi in the multivariate frailty models.}

\item{global_chisq}{Recurrent: a vector with the values of each multivariate Wald test.}
\item{dof_chisq}{Recurrent: a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test}{Recurrent: a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq}{Recurrent: a vector with the p_values for each global multivariate Wald test.}
\item{names.factor}{Recurrent: Names of the "as.factor" variables.}

\item{global_chisq_d}{Death: a vector with the values of each multivariate Wald test.}
\item{dof_chisq_d}{Death: a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test_d}{Death: a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq_d}{Death: a vector with the p_values for each global multivariate Wald test.}
\item{names.factordc}{Death: Names of the "as.factor" variables.}

\item{global_chisq2}{Recurrent2: a vector with the values of each multivariate Wald test.}
\item{dof_chisq2}{Recurrent2: a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test2}{Recurrent2: a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq2}{Recurrent2: a vector with the p_values for each global multivariate Wald test.}
\item{names.factor2}{Recurrent2: Names of the "as.factor" variables.}
}

\note{
"Kappa" (kappa[1], kappa[2] and kappa[3]) and "n.knots" (n.knots[1], n.knots[2] and n.knots[3]) are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}

\references{

Mazroui Y, Mathoulin-Pellissier S, MacGrogan G, Brouste V, Rondeau V. Multivariate frailty models for two types of recurrent events with an informative terminal event : Application to breast cancer data.
\emph{Biometrical journal}. In revision, 2013.

}
\author{
Virginie Rondeau, Yassin Mazroui, Alexandre Laurent and Amadou Diakite
}

\seealso{
 \code{\link{terminal}},\code{\link{event2}},
 \code{\link{print.multivePenal}},\code{\link{summary.multivePenal}},\code{\link{plot.multivePenal}}
}

\examples{
\dontrun{

data(dataMultiv)
## Calendar-time ##

MultivSpli <- multivePenal(Surv(t0,t1,deltaEvent1)~cluster(id)+var1+var2+
terminal(deltadc)+event2(deltaEvent2),formula.terminalEvent=~var1,
formula2=~var1+var2+var3,data=dataMultiv,Frailty=TRUE,recurrentAG=TRUE,
cross.validation=F,n.knots=c(4,4,4),kappa=c(1,1,1),hazard="Splines")

## print a fit
MultivSpli

## summary a fit
summary(MultivSpli)

## plot a fit
plot(MultivSpli,type.plot="haz",event="recurrent",conf.bands=TRUE)
}
}

\keyword{models}
\keyword{methods}
\keyword{multive}