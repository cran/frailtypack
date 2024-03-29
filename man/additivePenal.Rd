% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additivePenal.R
\name{additivePenal}
\alias{additivePenal}
\title{Fit an Additive Frailty model using a semiparametric penalized likelihood
estimation or a parametric estimation}
\usage{
additivePenal(formula, data, correlation = FALSE, recurrentAG =
FALSE, cross.validation = FALSE, n.knots, kappa, maxit = 350, hazard =
"Splines", nb.int, LIMparam = 1e-4, LIMlogl = 1e-4, LIMderiv = 1e-3,
print.times = TRUE)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a
\eqn{\sim} operator, and the terms on the right.  The response must be a
survival object as returned by the 'Surv' function like in survival package.
The \code{slope()} function is required.  Interactions are possible using *
or :.}

\item{data}{a 'data.frame' with the variables used in 'formula'.}

\item{correlation}{Logical value. Are the random effects correlated? If so,
the correlation coefficient is estimated. The default is FALSE.}

\item{recurrentAG}{Always FALSE for additive models (left-truncated data are
not allowed).}

\item{cross.validation}{Logical value. Is cross validation procedure used
for estimating smoothing parameter in the penalized likelihood estimation?
If so a search of the smoothing parameter using cross validation is done,
with kappa as the seed.  The cross validation is not implemented for two
strata. The default is FALSE.}

\item{n.knots}{integer giving the number of knots to use. Value required in
the penalized likelihood estimation.  It corresponds to the (n.knots+2)
splines functions for the approximation of the hazard or the survival
functions.  Number of knots must be between 4 and 20. (See Note)}

\item{kappa}{positive smoothing parameter in the penalized likelihood
estimation. In a stratified additive model, this argument must be a vector
with kappas for both strata. The coefficient kappa of the integral of the
squared second derivative of hazard function in the fit. To obtain an
initial value for \code{kappa}, a solution is to fit the corresponding
shared frailty model using cross validation (see cross.validation). We
advise the user to identify several possible tuning parameters, note their
defaults and look at the sensitivity of the results to varying them. Value
required. (See Note)}

\item{maxit}{maximum number of iterations for the Marquardtt algorithm.
Default is 350}

\item{hazard}{Type of hazard functions: "Splines" for semiparametric hazard
functions with the penalized likelihood estimation, "Piecewise-per" for
piecewise constant hazards functions using percentile, "Piecewise-equi" for
piecewise constant hazard functions using equidistant intervals, "Weibull"
for parametric Weibull functions. Default is "Splines".}

\item{nb.int}{Number of intervals (between 1 and 20) for the parametric
hazard functions ("Piecewise-per", "Piecewise-equi").}

\item{LIMparam}{Convergence threshold of the Marquardt algorithm for the
parameters (see Details), \eqn{10^{-4}} by default.}

\item{LIMlogl}{Convergence threshold of the Marquardt algorithm for the
log-likelihood (see Details), \eqn{10^{-4}} by default.}

\item{LIMderiv}{Convergence threshold of the Marquardt algorithm for the
gradient (see Details), \eqn{10^{-3}} by default.}

\item{print.times}{a logical parameter to print iteration process. Default
is TRUE.}
}
\value{
An additive model or more generally an object of class
'additivePenal'.  Methods defined for 'additivePenal' objects are provided
for print, plot and summary.

\item{b}{sequence of the corresponding estimation of the splines
coefficients, the random effects variances and the regression coefficients.}
\item{call}{The code used for fitting the model.} \item{coef}{the regression
coefficients.} \item{cov}{covariance between the two frailty terms
\eqn{(\bold{cov}(u_i,v_i))}} \item{cross.Val}{Logical value. Is cross
validation procedure used for estimating the smoothing parameters in the
penalized likelihood estimation?} \item{correlation}{Logical value. Are the
random effects correlated?} \item{DoF}{degrees of freedom associated with
the "kappa".}

\item{formula}{the formula part of the code used for the model.}
\item{groups}{the maximum number of groups used in the fit.} \item{kappa}{ A
vector with the smoothing parameters in the penalized likelihood estimation
corresponding to each baseline function as components.}
\item{loglikPenal}{the complete marginal penalized log-likelihood in the
semiparametric case.} \item{loglik}{the marginal log-likelihood in the
parametric case.} \item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit.}
\item{n.iter}{number of iterations needed to converge.} \item{n.knots
}{number of knots for estimating the baseline functions.}
\item{n.strat}{number of stratum.} \item{rho}{the corresponding correlation
coefficient for the two frailty terms.} \item{sigma2}{Variance for the
random intercept (the random effect associated to the baseline hazard
functions).} \item{tau2}{Variance for the random slope (the random effect
associated to the treatment effect across trials).} \item{varH}{the variance
matrix of all parameters before positivity constraint transformation
(Sigma2, Tau2, the regression coefficients and the spline coefficients).
Then after, the delta method is needed to obtain the estimated variance
parameters.} \item{varHIH}{the robust estimation of the variance matrix of
all parameters (Sigma2, Tau2, the regression coefficients and the spline
coefficients).} \item{varSigma2}{ The variance of the estimates of
"sigma2".} \item{varTau2}{ The variance of the estimates of "tau2".}
\item{varcov}{ Variance of the estimates of "cov".} \item{x}{matrix of times
where both survival and hazard functions are estimated. By default
seq(0,max(time),length=99), where time is the vector of survival times.}
\item{lam}{array (dim=3) of hazard estimates and confidence bands.}
\item{surv}{array (dim=3) of baseline survival estimates and confidence
bands.} \item{median}{The value of the median survival and its confidence bands. If there are
two stratas or more, the first value corresponds to the value for the 
first strata, etc.} \item{type.of.hazard}{Type of hazard functions (0:"Splines",
"1:Piecewise", "2:Weibull").} \item{type.of.Piecewise}{Type of Piecewise
hazard functions (1:"percentile", 0:"equidistant").} \item{nbintervR}{Number
of intervals (between 1 and 20) for the parametric hazard functions
("Piecewise-per", "Piecewise-equi").} \item{npar}{number of parameters.}
\item{nvar}{number of explanatory variables.} \item{noVar}{indicator of
explanatory variable.} \item{LCV}{the approximated likelihood
cross-validation criterion in the semiparametric case (with H minus the
converged Hessian matrix, and l(.) the full
log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{the Akaike information Criterion for the parametric
case.\deqn{AIC=\frac{1}{n}(np - l(.))}} \item{n.knots.temp}{initial value
for the number of knots.} \item{shape.weib}{shape parameter for the Weibull
hazard function.} \item{scale.weib}{scale parameter for the Weibull hazard
function.} \item{martingale.res}{martingale residuals for each cluster.}
\item{frailty.pred}{empirical Bayes prediction of the first frailty term.}
\item{frailty.pred2}{empirical Bayes prediction of the second frailty term.}
\item{linear.pred}{linear predictor: uses simply "Beta'X + u_i + v_i * X_1"
in the additive Frailty models.} \item{global_chisq}{a vector with the
values of each multivariate Wald test.} \item{dof_chisq}{a vector with the
degree of freedom for each multivariate Wald test.}
\item{global_chisq.test}{a binary variable equals to 0 when no multivariate
Wald is given, 1 otherwise.} \item{p.global_chisq}{a vector with the
p_values for each global multivariate Wald test.} \item{names.factor}{Names
of the "as.factor" variables.} \item{Xlevels}{vector of the values that
factor might have taken.} \item{contrasts}{type of contrast for factor
variable.} \item{beta_p.value}{p-values of the Wald test for the estimated
regression coefficients.}
}
\description{
{
Fit an additive frailty model using a semiparametric penalized likelihood
estimation or a parametric estimation. The main issue in a meta-analysis
study is how to take into account the heterogeneity between trials and
between the treatment effects across trials. Additive models are
proportional hazard model with two correlated random trial effects that act
either multiplicatively on the hazard function or in interaction with the
treatment, which allows studying for instance meta-analysis or multicentric
datasets.  Right-censored data are allowed, but not the left-truncated data.
A stratified analysis is possible (maximum number of strata = 2). This
approach is different from the shared frailty models.

\ifelse{html}{In an additive model, the hazard function for the j\out{<sup>th</sup>}
subject in the i\out{<sup>th</sup>} trial with random trial effect u\out{<sub>i</sub>} as
well as the random treatment-by-trial interaction v\out{<sub>i</sub>} is:

{\figure{additivemodel.png}{options: width="60\%"}}

where \eqn{\lambda}\out{<sub>0</sub>(0)} is the baseline hazard function, \eqn{\beta}\out{<sub>k</sub>} the
fixed effect associated to the covariate X\out{<sub>ijk</sub>} (k=1,..,p),
\eqn{\beta}\out{<sub>1</sub>} is the treatment effect and X\out{<sub>ij1</sub>} the treatment
variable. \eqn{\rho} is the corresponding correlation coefficient for the two frailty terms.}{In an additive model, the hazard function for the \eqn{j^{th}} subject in
the \eqn{i^{th}} trial with random trial effect \eqn{u_i} as well as the
random treatment-by-trial interaction \eqn{v_i} is:

\deqn{\left\{ \begin{array}{ll}
\lambda_{ij}(t|u_i,v_i)=\lambda_0(t)exp(u_i+v_iX_{ij1}+\sum_{k=1}^{p}\beta_kX_{ijk})
\\ \bold{cov}(u_i,v_i)=\bold{\rho\sigma\tau} \\
u_i\sim\bold{\mathcal{N}}(0,\bold{\sigma^2}) \mbox{,}
v_i\sim\bold{\mathcal{N}}(0,\bold{\tau^2}) \end{array} \right. }

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{\beta_k} the
fixed effect associated to the covariate \eqn{X_{ijk}} (k=1,..,p),
\eqn{\beta_1} is the treatment effect and \eqn{X_{ij1}} the treatment
variable. \eqn{\rho} is the corresponding correlation coefficient for the
two frailty terms.}
}
}
\details{
{
The estimated parameter are obtained by maximizing the penalized
log-likelihood or by a simple log-likelihood (in the parametric case) using
the robust Marquardtt algorithm (Marquardtt,1963).  The parameters are
initialized with values obtained with Cox proportional hazard model. The
iterations are stopped when the difference between two consecutive
loglikelhoods was small \eqn{(<10^{-4})}, the estimated coefficients were
stable (consecutive values \eqn{(<10^{-4})}, and the gradient small enough
\eqn{(<10^{-3})}. To be sure of having a positive function at all stages of
the algorithm, the spline coefficients were reparametrized to be positive at
each stage. The variance space of the two random effects is reduced, so the
variances are positive, and the correlation coefficient values are
constrained to be between -1 and 1.  The marginal log-likelihood depends on
integrations that are approximated by using the Laplace integration
technique with a first order approximation.  The smoothing parameter can be
fixed or estimated by maximizing likelihood cross-validation criterion. The
usual squared Wald statistic was modified to a mixture of two \eqn{\chi^2}
distribution to get significance test for the variance of the random
effects.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.1. An
adjusted Cox model is fitted, it provides new initial values for the splines
coefficients and the regression coefficients. The variances of the frailties
are initialized to 0.1. Then an additive frailty model with independent
frailties is fitted. At last, an additive frailty model with correlated
frailties is fitted.
}
}
\note{
"kappa" and "n.knots" are the arguments that the user have to change if the
fitted model does not converge.  "n.knots" takes integer values between 4
and 20. But with n.knots=20, the model would take a long time to converge.
So, usually, begin first with n.knots=7, and increase it step by step until
it converges. "kappa" only takes positive values. So, choose a value for
kappa (for instance 10000), and if it does not converge, multiply or divide
this value by 10 or 5 until it converges.
}
\examples{


\donttest{

###--- Additive model with 1 covariate ---###

data(dataAdditive)

modAdd <- additivePenal(Surv(t1,t2,event)~cluster(group)+
var1+slope(var1),correlation=TRUE,data=dataAdditive,
n.knots=8,kappa=10000)

#-- Var1 is boolean as a treatment variable

}


}
\references{
V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack:
An R package for the analysis of correlated survival data with frailty
models using penalized likelihood estimation or parametric estimation.
\emph{Journal of Statistical Software} \bold{47}, 1-28.

V. Rondeau, S. Michiels, B. Liquet, and J. P. Pignon (2008). Investigating
trial and treatment heterogeneity in an individual patient data
meta-analysis of survival data by mean of the maximum penalized likelihood
approach. \emph{Statistics in Medecine}, \bold{27}, 1894-1910.
}
\seealso{
\code{\link{slope}}
}
\keyword{file}
