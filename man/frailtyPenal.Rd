
\name{frailtyPenal}
\alias{frailtyPenal}
\alias{waldtest}
\alias{factor.names}
\alias{timedep.names}

\title{ Fit a Shared, Joint or Nested Frailty model }

\description{

\bold{Shared Frailty model}

Fit a shared gamma or log-normal frailty model using a semiparametric Penalized Likelihood estimation or parametric estimation on the hazard function.
Left truncated, right-censored data, interval-censored data and strata (max=2) are allowed.
It allows to obtain a non-parametric smooth hazard of survival function.
This approach is different from the partial penalized likelihood approach of Therneau et al.

The hazard function, conditionnal on the frailty term \eqn{\omega_i}, of a shared gamma frailty model for the \eqn{j^{th}} subject in the \eqn{i^{th}} group:

\deqn{\lambda_{ij}(t|\omega_i)=\lambda_0(t)\omega_i\exp(\bold{\beta^{'}Z_{ij}})}

\deqn{\omega_i\sim\Gamma\left(\frac{1}{\theta},\frac{1}{\theta}\right)
 \hspace{0.5cm} \bold{E}(\omega_i)=1 \hspace{0.5cm}\bold{Var}(\omega_i)=\theta}

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{\bold{\beta}} the vector of the regression coefficient associated to the covariate vector \eqn{\bold{Z_{ij}}} for the \eqn{j^{th}} individual in the \eqn{i^{th}} group.

Otherwise, in case of a shared log-normal frailty model, we have for the \eqn{j^{th}} subject in the \eqn{i^{th}} group:

\deqn{\lambda_{ij}(t|\eta_i)=\lambda_0(t)\exp(\eta_i+\bold{\beta^{'}Z_{ij}})}

\deqn{\eta_i\sim N(0,\sigma^2)}

From now on, you can also consider time-varying effects covariates in your model, see \code{timedep} function for more details.

\bold{Joint Frailty model}

Fit a joint either with gamma or log-normal frailty model for recurrent and terminal events using a penalized likelihood estimation on the hazard function or a parametric estimation. Right-censored data is allowed. Left-truncated data and stratified analysis are not possible. Joint frailty models allow studying, jointly, survival processes of recurrent and terminal events, by considering the terminal event as an informative censoring. 

There is two kinds of joint frailty models that can be fitted with \code{frailtyPenal} :

- The first one (Rondeau et al. 2007) includes a common frailty term to the individuals \eqn{(\omega_i)} for the two rates which will take into account the heterogeneity in the data, associated with unobserved covariates. The frailty term acts differently for the two rates ( \eqn{\omega_i} for the recurrent rate and \eqn{\omega_i^{\alpha}} for the death rate). The covariates could be different for the recurrent rate and death rate.

For the \eqn{j^{th}}{j^th} recurrence \eqn{(j=1,...,n_i)} and the \eqn{i^{th}}{i^th} subject \eqn{(i=1,...,G)}, the joint gamma frailty model for recurrent event hazard function \eqn{r_{ij}(.)} and death rate \eqn{\lambda_i(.)} is :

\deqn{\left\{
\begin{array}{ll}
r_{ij}(t|\omega_i)=\omega_ir_0(t)\exp(\bold{\beta_1^{'}Z_i(t)}) & \mbox{(Recurrent)} \\
\lambda_i(t|\omega_i)=\omega_i^{\alpha}\lambda_0(t)\exp(\bold{\beta_2^{'}Z_i(t)}) & \mbox{(Death)} \\
 \end{array}
 \right.
}

where \eqn{r_0(t)} (resp. \eqn{\lambda_0(t)}) is the recurrent (resp. terminal) event baseline hazard function, \eqn{\bold{\beta_1}} (resp. \eqn{\bold{\beta_2}}) the regression coefficient vector, \eqn{\bold{Z_i(t)}} the covariate vector. The random effects of frailties \eqn{\omega_i\sim\bold{\Gamma}(\frac{1}{\theta},\frac{1}{\theta})} and are iid.

The joint log-normal frailty model will be :

\deqn{\left\{
\begin{array}{ll}
r_{ij}(t|\eta_i)=r_0(t)\exp(\eta_i+\bold{\beta_1^{'}Z_i(t)}) & \mbox{(Recurrent)} \\
\lambda_i(t|\eta_i)=\lambda_0(t)\exp(\alpha \eta_i+\bold{\beta_2^{'}Z_i(t)}) & \mbox{(Death)} \\
 \end{array}
 \right.
}

where \deqn{\eta_i\sim N(0,\sigma^2)}

- The second one (Rondeau et al. 2011) is quite similar but the frailty term is common to the individuals from a same group. This model is useful for joint modelling two clustered survival outcomes. In this case, j is for the subject and i for the cluster.

\deqn{\left\{
\begin{array}{ll}
r_{ij}(t|u_i)=u_ir_0(t)\exp(\bold{\beta_1^{'}Z_{ij}(t)}) & \mbox{(Time to event)} \\
\lambda_{ij}(t|u_i)=u_i^{\alpha}\lambda_0(t)\exp(\bold{\beta_2^{'}Z_{ij}(t)}) & \mbox{(Death)} \\
 \end{array}
 \right.
}

In case of a log-normal distribution of the frailties, we will have :

\deqn{\left\{
\begin{array}{ll}
r_{ij}(t|v_i)=r_0(t)\exp(v_i+\bold{\beta_1^{'}Z_{ij}(t)}) & \mbox{(Time to event)} \\
\lambda_{ij}(t|v_i)=\lambda_0(t)\exp(\alpha v_i+\bold{\beta_2^{'}Z_{ij}(t)}) & \mbox{(Death)} \\
 \end{array}
 \right.
}

where \deqn{v_i\sim N(0,\sigma^2)}

This joint frailty model can also be applied to clustered recurrent events and a terminal event (example on "readmission" data below).

From now on, you can also consider time-varying effects covariates in your model, see \code{timedep} function for more details.

\bold{Nested Frailty model}

Fit a nested frailty model using a Penalized Likelihood on the hazard function or using a parametric estimation. Nested frailty models allow survival studies for hierarchically clustered data by including two iid gamma random effects.
Left truncated and censored data are allowed. Stratification analysis is allowed(maximum of strata=2).

The hazard function conditional on the two frailties \eqn{v_i} and \eqn{w_{ij}} for the \eqn{k^{th}} individual of the \eqn{j^{th}} subgroup of the \eqn{i^{th}} group is :

\deqn{\left\{
\begin{array}{ll}
\lambda_{ijk}(t|v_i,w_{ij})=v_iw_{ij}\lambda_0(t)exp(\bold{\beta^{'}X_{ijk}}) \\
v_i\sim\Gamma\left(\frac{1}{\alpha},\frac{1}{\alpha}\right) \hspace{0.05cm}i.i.d. \hspace{0.2cm} \bold{E}(v_i)=1 \hspace{0.2cm}\bold{Var}(v_i)=\alpha \\
 w_{ij}\sim\Gamma\left(\frac{1}{\eta},\frac{1}{\eta}\right)\hspace{0.05cm}i.i.d. \hspace{0.2cm} \bold{E}(w_{ij})=1 \hspace{0.2cm} \bold{Var}(w_{ij})=\eta
 \end{array}
 \right.
}

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{X_{ijk}}
 denotes the covariate vector and \eqn{\beta} the corresponding vector of regression parameters.

}

\usage{

frailtyPenal(formula, formula.terminalEvent, data, 
             recurrentAG = FALSE, cross.validation = FALSE,
             n.knots, kappa1, kappa2, maxit = 350, hazard =
             "Splines", nb.int1, nb.int2, RandDist = "Gamma",
             betaknots = 1, betaorder = 3, init.B, init.Theta,
             init.Alpha, Alpha, LIMparam = 1e-4, LIMlogl = 1e-4,
             LIMderiv = 1e-4, print.times = TRUE)
}

\arguments{

  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right. The response must be a
          survival object as returned by the 'Surv' function like in survival package.
          In case of interval-censored data, the response must be an object as returned by
          the 'SurvIC' function from this package }

  \item{formula.terminalEvent}{only for joint model: a formula object, only requires terms on the right to indicate which variables are modelling the terminal event.}

  \item{data}{a 'data.frame' with the variables used in 'formula'.}

  \item{recurrentAG}{Logical value. Is Andersen-Gill model fitted? If so indicates that recurrent event times with the
                 counting process approach of Andersen and Gill is used. This formulation can be used for dealing with
                 time-dependent covariates. The default is FALSE.}

  \item{cross.validation}{Logical value. Is cross validation procedure used for estimating smoothing parameter
                 in the penalized likelihood estimation?
                 If so a search of the smoothing parameter using cross validation is done, with kappa1 as the seed.
                 The cross validation is not implemented for two strata, neither for interval-censored data. The cross validation has been implemented 
                 for a Cox proportional hazard model, with no covariates. The default is FALSE.}

  \item{n.knots}{ integer giving the number of knots to use. Value required in the penalized likelihood estimation.
                It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions.
                We estimate I or M-splines of order 4. When the user set a number of knots equals to k (n.knots=k) then the number of interior knots is (k-2) and the number of splines is (k-2)+order.  
                Number of knots must be between 4 and 20. (See Note)}

  \item{kappa1}{ positive smoothing parameter in the penalized likelihood estimation. The coefficient kappa of the integral of the squared
          second derivative of hazard function in the fit (penalized log likelihood). To obtain an initial 
          value for \code{kappa1} (or \code{kappa2}), a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation).
          We advise the user to identify several possible tuning parameters, note their defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

  \item{kappa2}{positive smoothing parameter in the penalized likelihood estimation for the second stratum when data are stratified or for the terminal event rate for joint model (stratification not allowed in joint model). See kappa1. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

  \item{hazard}{Type of hazard functions: "Splines" for semiparametric hazard functions with the penalized likelihood estimation, 
            "Piecewise-per" for piecewise constant hazard function using percentile (not available for interval-censored data), "Piecewise-equi" for piecewise constant hazard 
            function using equidistant intervals, "Weibull" for parametric weibull functions. Default is "Splines".}

  \item{nb.int1}{Number of intervals (between 1 and 20) for the parametric hazard functions ("Piecewise-per", "Piecewise-equi"). In a joint model, this argument is used for the recurrent hazard function.}

  \item{nb.int2}{Only for joint model : number of intervals (between 1 and 20) for the death parametric hazard function ("Piecewise-per", "Piecewise-equi").}

  \item{RandDist}{Type of random effect distribution: "Gamma" for a gamma distribution, "LogN" for a log-normal distribution. Default is "Gamma". Not implemented for nested model.}

  \item{betaknots}{Number of inner knots used for the estimation of B-splines. Default is 1. See 'timedep' function for more details. Not implemented for nested model.}

  \item{betaorder}{Order of the B-splines. Defaut is cubic B-splines (order = 3). See 'timedep' function for more details. Not implemented for nested model.}
  
  \item{init.B}{A vector of inital values for regression coefficients. This vector should be of the same size as the whole vector of covariates. Default is 0.1 for each (for Cox and shared model) or 0.5 (for joint model).}
  
  \item{init.Theta}{Inital value for variance of the frailties.}
  
  \item{init.Alpha}{Only for joint model : initial value for parameter alpha.}

  \item{Alpha}{Only for joint model : input "none" so as to fit a joint model without the parameter alpha.}

  \item{LIMparam}{Convergence threshold of the Marquard algorithm for the parameters (see Details), \eqn{10^{-4}} by default.}
  
  \item{LIMlogl}{Convergence threshold of the Marquard algorithm for the log-likelihood (see Details), \eqn{10^{-4}} by default.}
  
  \item{LIMderiv}{Convergence threshold of the Marquard algorithm for the gradient (see Details), \eqn{10^{-4}} by default.}

  \item{print.times}{a logical parameter to print iteration process. Default is TRUE.}
}

\details{
Typical usages are
for a Cox model
\preformatted{frailtyPenal(Surv(time,event)~var1+var2, data, \dots)}

for a shared model
\preformatted{frailtyPenal(Surv(time,event)~cluster(group)+var1+var2,
             data, \dots)}

for a joint model
\preformatted{frailtyPenal(Surv(time,event)~cluster(group)+var1+var2+
             var3+terminal(death), formula.terminalEvent=~
             var1+var4, data, \dots)}

for a joint model for clustered data
\preformatted{frailtyPenal(Surv(time,event)~cluster(group)+num.id(group2)+
             var1+var2+var3+terminal(death),
             formula.terminalEvent=~var1+var4, data, \dots)}

for a nested model
\preformatted{frailtyPenal(Surv(time,event)~cluster(group)+subcluster(sbgroup)+
             var1+var2, data, \dots)}

The estimated parameter are obtained using the robust Marquardt algorithm (Marquardt, 1963) which is a combination
between a Newton-Raphson algorithm and a steepest descent algorithm. The iterations are stopped when the difference
between two consecutive loglikelihoods was small \eqn{(<10^{-4})}, the estimated coefficients were stable
(consecutive values \eqn{(<10^{-4})}, and the gradient small enough \eqn{(<10^{-4})}. When frailty parameter
is small, numerical problems may arise. To solve this problem, an alternative formula of the penalized
log-likelihood is used (see Rondeau, 2003 for further details). Cubic M-splines of order 4 are used for
the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.

The inverse of the hessian matrix is the variance estimator and to deal with the positivity constraint of
the variance component and the spline coefficients, a squared transformation is used and the standard errors
are computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth parameter can be chosen by
maximizing a likelihood cross validation criterion (Joly and other, 1998). The integrations in the full
log likelihood were evaluated using Gaussian quadrature. Laguerre polynomials with 20 points were used
to treat the integrations on \eqn{[0,\infty[}

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.1. In case of shared model, the program fits,
firstly, an adjusted Cox model to give new initial values for the splines and the regression coefficients.
The variance of the frailty term \eqn{\theta} is initialized to 0.1. Then, a shared frailty model is fitted.

In case of a joint frailty model, the splines and the regression coefficients are initialized to 0.5. The
program fits an adjusted Cox model to have new initial values for the regression and the splines coefficients.
The variance of the frailty term \eqn{\theta} and the coefficient \eqn{\alpha} associated in the death hazard
function are initialized to 1. Then, it fits a joint frailty model.

In case of a nested model, the program fits an adjusted Cox model to provide new initial values for the
regression and the splines coefficients. The variances of the frailties are initialized to 0.1. Then,
a shared frailty model with covariates with only subgroup frailty is fitted to give a  new initial value
for the variance of the subgroup frailty term. Then, a shared frailty model with covariates and only group
frailty terms is fitted to give a  new initial value for the variance of the group frailties. In a last step,
a nested frailty model is fitted.

}

\value{

The following components are included in a 'frailtyPenal' object for each model.

\item{b}{sequence of the corresponding estimation of the coefficients for the hazard functions (parametric or semiparametric), the random effects variances and the regression coefficients.}
\item{call}{The code used for the model.}
\item{formula}{the formula part of the code used for the model.}
\item{coef}{the regression coefficients.}
\item{cross.Val }{Logical value. Is cross validation procedure used for estimating the smoothing parameters in the penalized likelihood estimation?}
\item{DoF}{Degrees of freedom associated with the "kappa".}
\item{formula }{the formula part of the code used for the model.}
\item{groups}{the maximum number of groups used in the fit.}
\item{kappa}{ A vector with the smoothing parameters  in the penalized likelihood estimation corresponding to each baseline function as components.}
\item{lam}{matrix of hazard estimates and confidence bands.}
\item{lam2}{the same value as lam for the second stratum.}
\item{loglikPenal}{the complete marginal penalized log-likelihood in the semiparametric case.}
\item{loglik}{the marginal log-likelihood in the parametric case.}
\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit.}
\item{n.iter}{number of iterations needed to converge.}
\item{n.knots}{number of knots for estimating the baseline functions in the penalized likelihood estimation.}
\item{n.strat}{ number of stratum.}
\item{surv}{matrix of baseline survival estimates and confidence bands.}
\item{surv2}{the same value as surv for the second stratum.}
\item{varH}{the variance matrix of all parameters before positivity constraint transformation. Thenafter, the delta method is needed to obtain the estimated variance parameters. That is why some variances don't match with the printed values at the end of the model.}
\item{varHIH}{the robust estimation of the variance matrix of all parameters.}
\item{x1}{vector of times where both survival and hazard function are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{vector of times for the second stratum or the terminal event (see x1 value).}
\item{type.of.hazard}{Type of hazard functions (0:"Splines", "1:Piecewise", "2:Weibull").}
\item{type.of.Piecewise}{Type of Piecewise hazard functions (1:"percentile", 0:"equidistant").}
\item{nbintervR}{Number of intervals (between 1 and 20) for the parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{npar}{number of parameters.}
\item{nvar}{number of explanatory variables.}
\item{noVar}{indicator of explanatory variables.}
\item{LCV}{the approximated likelihood cross-validation criterion in the semiparametric case (with H minus the converged hessien matrix, and l(.) the full log-likelihood).\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{the Akaike information Criterion for the parametric case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
\item{n.knots.temp}{initial value for the number of knots.}
\item{shape.weib}{shape parameter for the Weibull hazard function.}
\item{scale.weib}{scale parameter for the Weibull hazard function.}
\item{martingale.res}{martingale residuals for each cluster.}
\item{martingaleCox}{martingale residuals for observation in the Cox model.}
\item{Frailty}{Logical value. Was model with frailties fitted ?}
\item{frailty.pred}{empirical Bayes prediction of the frailty term (ie, using conditional posterior distributions).}
\item{frailty.var}{variance of the empirical Bayes prediction of the frailty term (only for gamma frailty models).}
\item{frailty.sd}{standard error of the frailty empirical Bayes prediction (only for gamma frailty models).}
\item{global_chisq}{a vector with the values of each multivariate Wald test.}
\item{dof_chisq}{a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test}{a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq}{a vector with the p_values for each global multivariate Wald test.}
\item{names.factor}{Names of the "as.factor" variables.}
\item{Xlevels}{vector of the values that factor might have taken.}
\item{contrasts}{type of contrast for factor variable.}

The following components are specific to \bold{shared} models.

\item{theta}{variance of the gamma frailty parameter \eqn{(\bold{Var}(\omega_i))}}
\item{sigma2}{variance of the log-normal frailty parameter \eqn{(\bold{Var}(\eta_i))}}
\item{linear.pred}{linear predictor: uses simply "Beta'X" in the cox proportional hazard model or "Beta'X + log w_i" in the shared gamma frailty models, otherwise uses "Beta'X + w_i" for log-normal frailty distribution.}
\item{BetaTpsMat}{matrix of time varying-effects and confidence bands (the first column used for absciss of times)}

The following components are specific to \bold{joint} models.

\item{theta}{variance of the gamma frailty parameter \eqn{(\bold{Var}(\omega_i))} or \eqn{(\bold{Var}(u_i))}}
\item{sigma2}{variance of the log-normal frailty parameter \eqn{(\bold{Var}(\eta_i))} or \eqn{(\bold{Var}(v_i))}}
\item{indic_alpha}{indicator if a joint frailty model with \eqn{\alpha} parameter was fitted}
\item{alpha}{the coefficient \eqn{\alpha} associated with the frailty parameter in the terminal hazard function.}
\item{nbintervR}{Number of intervals (between 1 and 20) for the recurrent parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nbintervDC}{Number of intervals (between 1 and 20) for the death parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nvar}{A vector with the number of covariates of each type of hazard function as components.}
\item{nvarRec}{number of recurrent explanatory variables.}
\item{nvarEnd}{number of death explanatory variables.}
\item{noVar1}{indicator of recurrent explanatory variables.}
\item{noVar2}{indicator of death explanatory variables.}
\item{martingale.res}{martingale residuals for each cluster (recurrent).}
\item{martingaledeath.res}{martingale residuals for each cluster (death).}
\item{linear.pred}{linear predictor: uses "Beta'X + log w_i" in the gamma frailty model, otherwise uses "Beta'X + eta_i" for log-normal frailty distribution}
\item{lineardeath.pred}{linear predictor for the terminal part : "Beta'X + alpha.log w_i" for gamma, "Beta'X + alpha.eta_i" for log-normal frailty distribution}
\item{Xlevels}{vector of the values that factor might have taken for the recurrent part.}
\item{contrasts}{type of contrast for factor variable for the recurrent part.}
\item{Xlevels2}{vector of the values that factor might have taken for the death part.}
\item{contrasts2}{type of contrast for factor variable for the death part.}
\item{BetaTpsMat}{matrix of time varying-effects and confidence bands for recurrent event (the first column used for absciss of times of recurrence)}
\item{BetaTpsMatDc}{matrix of time varying-effects and confidence bands for terminal event (the first column used for absciss of times of death)}

The following components are specific to \bold{nested} models.

\item{alpha}{variance of the cluster effect \eqn{(\bold{Var}(v_{i}))}}
\item{eta}{variance of the subcluster effect \eqn{(\bold{Var}(w_{ij}))}}
\item{subgroups}{the maximum number of subgroups used in the fit.}
\item{frailty.pred.group}{empirical Bayes prediction of the frailty term by group.}
\item{frailty.pred.subgroup}{empirical Bayes prediction of the frailty term by subgroup.}
\item{linear.pred}{linear predictor: uses "Beta'X + log v_i.w_ij".}
\item{subgbyg}{subgroup by group.}
\item{n.strat}{A vector with the number of covariates of each type of hazard function as components.}

}

\note{

From a prediction aim, we recommand you to input a data sorted by the group variable with numerical numbers
from 1 to n (number of groups). In case of a nested model, we recommand you to input a data sorted by the
group variable then sorted by the subgroup variable both with numerical numbers from 1 to n (number of groups)
and from 1 to m (number or subgroups).
"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model
does not converge.
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to
converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does
not converge, multiply or divide this value by 10 or 5 until it converges.

}

\references{

V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack: An R package for the analysis of correlated survival data with frailty models using penalized likelihood estimation or parametric estimation.
\emph{Journal of Statistical Software} \bold{47}, 1-28.

V. Rondeau, J.P. Pignon, S. Michiels (2011). A joint model for the dependance between clustered times to tumour progression and deaths: A meta-analysis of chemotherapy in head and neck cancer.
\emph{Statistical methods in medical research} \bold{897}, 1-19.

V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P. Soubeyran (2007). Joint frailty models for recurring events and death using maximum penalized likelihood estimation:application on cancer events.
\emph{Biostatistics} \bold{8},4, 708-721.

V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum penalized likelihood estimation.
\emph{Statistics in Medecine}, \bold{25}, 4036-4052.

V. Rondeau, D. Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
\emph{Lifetime Data Analysis} \bold{9}, 139-153.

McGilchrist CA and Aisbett CW (1991). Regression with frailty in survival analysis.
\emph{Biometrics} \bold{47}, 461-466.

D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters.
\emph{SIAM Journal of Applied Mathematics}, 431-441.

}

\seealso{
 \code{\link{SurvIC}}, \code{\link{subcluster}}, \code{\link{terminal}}, \code{\link{num.id}}, \code{\link{timedep}}
}


\examples{

\dontrun{

###---  COX proportional hazard model (SHARED without frailties) ---###
###---  estimated with penalized likelihood ---###

data(kidney)
frailtyPenal(Surv(time,status)~sex+age,
n.knots=12,kappa1=10000,data=kidney)

###---  Shared Frailty model  ---###

frailtyPenal(Surv(time,status)~cluster(id)+sex+age,
n.knots=12,kappa1=10000,data=kidney)

#-- with an initialisation of regression coefficients

frailtyPenal(Surv(time,status)~cluster(id)+sex+age,
n.knots=12,kappa1=10000,data=kidney,init.B=c(-1.44,0))

#-- with truncated data

data(dataNested)

frailtyPenal(Surv(t1,t2,event) ~ cluster(group),
data=dataNested,n.knots=10,kappa1=10000,
cross.validation=TRUE,recurrentAG=FALSE)

#-- stratified analysis

data(readmission)
frailtyPenal(Surv(time,event)~cluster(id)+dukes+strata(sex),
n.knots=10,kappa1=10000,kappa2=10000,data=readmission)

#-- recurrentAG=TRUE

frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
charlson,data=readmission,n.knots=6,kappa1=1e5,recurrentAG=TRUE)

#-- cross.validation=TRUE

frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
charlson,data=readmission,n.knots=6,kappa1=5000,recurrentAG=TRUE,
cross.validation=TRUE)

#-- log-normal distribution

frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex+dukes+
charlson,data=readmission,n.knots=6,kappa1=5000,recurrentAG=TRUE,
RandDist="LogN")

###--- Joint Frailty model (recurrent and terminal events) ---###

data(readmission)
#-- Gap-time
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
terminal(death),formula.terminalEvent=~sex+dukes+charlson,
data=readmission,n.knots=10,kappa1=100,kappa2=100,
recurrentAG=FALSE)

#-- Calendar time
modJoint.calendar <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
+dukes+charlson+terminal(death),formula.terminalEvent=~sex
+dukes+charlson,data=readmission,n.knots=10,kappa1=100,
kappa2=100,recurrentAG=TRUE)

#-- without alpha parameter
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
terminal(death),formula.terminalEvent=~sex+dukes+charlson,
data=readmission,n.knots=10,kappa1=100,kappa2=100,
recurrentAG=FALSE,Alpha="none")

#-- log-normal distribution

modJoint.log <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
+dukes+charlson+terminal(death),formula.terminalEvent=~sex
+dukes+charlson,data=readmission,n.knots=10,kappa1=100,
kappa2=100,recurrentAG=TRUE,RandDist="LogN")

###--- Joint Frailty model for clustered data ---###

#-- here is a generated cluster (31 clusters of 13 subjects)
readmission <- transform(readmission,group=id%%31+1)

joi.clus <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(group)+
num.id(id)+dukes+charlson+sex+chemo+terminal(death),
formula.terminalEvent=~dukes+charlson+sex+chemo,
data=readmission,recurrentAG=TRUE,n.knots=10,
kappa1=100,kappa2=100)

###--- Nested Frailty model ---###

data(dataNested)
modClu <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,data=dataNested,
n.knots=8,kappa1=50000)

modClu.str <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+strata(cov2),data=dataNested,
n.knots=8,kappa1=50000,kappa2=50000)

}

}

\keyword{models}

