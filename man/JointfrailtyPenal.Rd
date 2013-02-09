
\name{frailtyPenal for Joint frailty models}
\alias{frailtyPenal for Joint frailty models}
\alias{frailtyPenal for Joint models}


\title{ Fit Joint Frailty model for recurrent and terminal events using semiparametric penalized likelihood estimation or a parametric estimation }

\description{
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

}


\synopsis{

frailtyPenal(formula, formula.terminalEvent, data, Frailty = FALSE,
            joint = FALSE, recurrentAG = FALSE, cross.validation =
            FALSE, n.knots, kappa1, kappa2, maxit = 350,
            hazard = "Splines", nb.int1, nb.int2, RandDist = "Gamma")
}

\arguments{

  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow}
          operator, and the terms on the right. The response must be a
          survival object as returned by the 'Surv' function like in survival package.
          To fit a joint model for clustered data, the \code{num.id()} function is required.}

  \item{formula.terminalEvent}{a formula object, only requires terms on the right to indicate which variables are modelling the terminal event.}

  \item{data}{a 'data.frame' in which to interpret the variables named in the 'formula' and 'formula.terminalEvent'.}

  \item{Frailty}{Logical value. Is model with frailties fitted? If so, variance of frailty parameter is estimated. The default is FALSE}

  \item{joint}{Logical value. Is joint model fitted? If so, 'formula.terminalEvent' is required. The default is FALSE}

  \item{recurrentAG}{Logical value. Is Andersen-Gill model fitted? If so indicates that recurrent event times with the 
                 counting process approach of Andersen and Gill is used. This formulation can be used for dealing with
                 time-dependent covariates. The default is FALSE.}

  \item{cross.validation}{Logical value. Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation?
                 If so a search of the smoothing parameter using cross validation is done, with kappa1 as the seed.
                 The cross validation is not implemented for two strata. The cross validation has been implemented 
                 for a Cox proportional hazard model, with no covariates. The default is FALSE.}

  \item{n.knots}{integer giving the number of knots to use. Value required in the penalized likelihood estimation.
                It corresponds to the (n.knots+2) splines functions for the approximation of the hazard or the survival functions.
                Number of knots must be between 4 and 20. (See Note)}

  \item{kappa1}{positive smoothing parameter in the penalized likelihood estimation. The coefficient kappa of the integral of the squared
              second derivative of hazard function in the fit (penalized log likelihood). To obtain a good
              value for \code{kappa1}, a solution is to fit the corresponding shared frailty model using cross validation
              (See cross.validation) with the event indicator as the event of interest. We advise the user to identify several 
              possible tuning parameters, note their defaults and look at the sensitivity of the results to varying them. Value required.(See Note)}

  \item{kappa2}{positive smoothing parameter in the penalized likelihood estimation for the terminal event rate. To obtain a good value
              for \code{kappa2}, a solution is to fit the corresponding Cox model using cross-validation (See cross.validation) with the
              death indicator as the event of interest. See kappa1. Stratification is not allowed here. }

  \item{maxit}{maximum number of iterations for the Marquardt algorithm. Default is 350}

  \item{hazard}{Type of hazard functions: "Splines" for semiparametric hazard functions with the penalized likelihood estimation,
              "Piecewise-per" for piecewise constant hazard function using percentile, "Piecewise-equi" for piecewise constant 
               hazard function using equidistant intervals, "Weibull" for parametric weibull function. Default is "Splines".}

  \item{nb.int1}{Number of intervals (between 1 and 20) for the recurrent parametric hazard functions ("Piecewise-per", "Piecewise-equi").}

  \item{nb.int2}{Number of intervals (between 1 and 20) for the death parametric hazard functions ("Piecewise-per", "Piecewise-equi").}

  \item{RandDist}{Type of random effect distribution: "Gamma" for a gamma distribution, "LogN" for a log-normal distribution. Default is "Gamma".}

}


\details{

  The estimated parameter are obtained by maximizing the semiparametric penalized loglikelihood estimation or a parametric estimation using the robust Marquardt algorithm (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm and a steepest descent algorithm. When frailty parameter is small, numerical problems may arise. To solve this problem, an alternative formula of the penalized log-likelihood is used (see Rondeau, 2003 for further details). 
  Cubic M-splines of order 4 are used for the hazard function, and I-splines (integrated M-splines) are used for the cumulative hazard function.
The inverse of the hessian matrix is the variance estimator and to deal with thepositivity constraint of the variance component, a squared transformation is used and the Standard errors are computed by the \eqn{\Delta}-method (Knight & Xekalaki, 2000). The smooth parameter can be chosen by maximizing a likelihood cross validation criterion (Joly and other, 1998). The integrations in the full log likelihood were evaluated using Gaussian quadrature. Laguerre polynomials with 20 points were used to treat the integrations on \eqn{[0,\infty[}.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.5. The program fits an adjusted Cox model to have new initial values for the regression and the splines coefficients. The variance of the frailty term \eqn{\theta} and the coefficient \eqn{\alpha} associated in the death hazard function are initialized to 1. Then, it fits a joint frailty model.

}

\value{
Parameters estimates of a joint frailty model, more generally a 'fraityPenal' object. Methods defined for 'frailtyPenal' objects are
provided for print,  plot and summary. The following components are included in a 'frailtyPenal' object for Joint frailty models.

\item{b}{sequence of the corresponding estimation of the splines coefficients, the random effects variances, the power of the frailties and the regression coefficients.}

\item{alpha}{the coefficient associated with the frailty parameter in the terminal hazard function.}
\item{call}{The code used for fitting the model.}
\item{coef}{the regression coefficients.}

\item{cross.Val}{Logical value. Is cross validation procedure used for estimating the smoothing parameters  in the penalized likelihood estimation?}

\item{formula}{the formula part of the code used for the model.}
\item{groups}{the maximum number of groups used in the fit.}
\item{kappa}{ A vector with the smoothing parameters  in the penalized likelihood estimation corresponding to each baseline function as components.}
\item{lam}{matrix of hazard estimates and confidence bands.}
\item{lam2}{the same value as lam for  the second stratum.}
\item{loglikPenal}{the complete marginal penalized log-likelihood in the semiparametric case.}
\item{loglik}{the marginal log-likelihood in the parametric case.}
\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit.}
\item{n.knots}{number of knots for estimating the baseline functions.}
\item{n.iter}{number of iterations needed to converge.}
\item{surv}{matrix of baseline survival estimates and confidence bands.}
\item{surv2}{the same value as surv for the  the second stratum.}

\item{theta}{variance of the gamma frailty parameter \eqn{(\bold{Var}(\omega_i))} or \eqn{(\bold{Var}(u_i))}}
\item{sigma2}{variance of the log-normal frailty parameter \eqn{(\bold{Var}(\eta_i))} or \eqn{(\bold{Var}(v_i))}}
 \item{varH}{the variance matrix of all parameters before positivity constraint transformation (theta, the regression coefficients and the spline coefficients). Thenafter, the delta method is needed to obtain the estimated variance parameters.}
 \item{varHIH}{the robust estimation of the variance matrix of all parameters (theta, the regression coefficients and the spline coefficients).}
\item{x1}{vector of times where both survival and hazard functions for the recurrent events are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
\item{x2}{vector of times for the second stratum (see x1 value).}
\item{type.of.hazard}{Type of hazard functions (0:"Splines", "1:Piecewise", "2:Weibull").}
\item{type.of.Piecewise}{Type of Piecewise hazard functions (1:"percentile", 0:"equidistant").}
\item{nbintervR}{Number of intervals (between 1 and 20) for the recurrent parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{nbintervDC}{Number of intervals (between 1 and 20) for the death parametric hazard functions ("Piecewise-per", "Piecewise-equi").}
\item{npar}{number of parameters.}
\item{nvar}{A vector with the number of covariates of each type of hazard function as components.}
\item{nvarRec}{number of recurrent explanatory vairables.}
\item{nvarEnd}{number of death explanatory vairables.}
\item{noVar1}{indicator of recurrent explanatory vairables.}
\item{noVar2}{indicator of death explanatory vairables.}
\item{LCV}{the approximated likelihood cross-validation criterion in the semiparametric case (with H minus the converged hessien matrix, and l(.) the full log-likelihood.\deqn{LCV=\frac{1}{n}(trace(H^{-1}_{pl}H) - l(.))}}
\item{AIC}{the Akaike information Criterion for the parametric case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
\item{n.knots.temp}{initial value for the number of knots.}
\item{shape.weib}{shape parameters for the weibull hazard function.}
\item{scale.weib}{scale parameters for the weibull hazard function.}
\item{martingale.res}{martingale residuals for each cluster (recurrent).}
\item{martingaledeath.res}{martingale residuals for each cluster (death).}
\item{frailty.pred}{empirical Bayes prediction of the frailty term.}
\item{linear.pred}{linear predictor: uses "Beta'X + log w_i" in the joint gamma frailty models, otherwise uses "Beta'X + eta_i" in log-normal}
\item{lineardeath.pred}{linear predictor for the terminal part form the joint model: "Beta'X + alpha.log w_i" for gamma, "Beta'X + alpha.eta_i" .}

\item{global_chisq}{a vector with the values of each multivariate Wald test.}
\item{dof_chisq}{a vector with the degree of freedom for each multivariate Wald test.}
\item{global_chisq.test}{a binary variable equals to 0 when no multivariate Wald is given, 1 otherwise.}
\item{p.global_chisq}{a vector with the p_values for each global multivariate Wald test.}
\item{names.factor}{Names of the "as.factor" variables.}
\item{Xlevels}{vector of the values that factor might have taken for the recurrent part.}
\item{contrasts}{type of contrast for factor variable for the recurrent part.}
\item{Xlevels2}{vector of the values that factor might have taken for the death part.}
\item{contrasts2}{type of contrast for factor variable for the death part.}
}

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}


\references{
V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack: An R package for the analysis of correlated survival data with frailty models using penalized likelihood estimation or parametric estimation.
\emph{Journal of Statistical Software} \bold{47}, 1-28.

V. Rondeau, J.P. Pignon, S. Michiels (2011). A joint model for the dependance between clustered times to tumour progression and deaths: A meta-analysis of chemotherapy in head and neck cancer.
\emph{Statistical methods in medical research} \bold{897}, 1-19.

V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P. Soubeyran (2007). Joint frailty models for recurring events and death using maximum penalized likelihood estimation:application on cancer events.
\emph{Biostatistics} \bold{8},4, 708-721.

V. Rondeau, D. Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model.
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

\dontrun{

data(readmission)

### Joint Frailty model ###

## Gap-time ##
modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
                  terminal(death),formula.terminalEvent=~sex+dukes+charlson,
                  data=readmission,n.knots=14,kappa1=9.55e+9,kappa2=1.41e+12,
                  Frailty=TRUE,joint=TRUE,recurrentAG=FALSE)

## Calendar time ##
modJoint.calendar <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
              +dukes+charlson+terminal(death),formula.terminalEvent=~sex
              +dukes+charlson,data=readmission,n.knots=10,kappa1=9.55e9,
              kappa2=1.41e12,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE)

### Joint Frailty model with log-normal distribution ###

modJoint.log <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+sex
              +dukes+charlson+terminal(death),formula.terminalEvent=~sex
              +dukes+charlson,data=readmission,n.knots=10,kappa1=9.55e9,
              kappa2=1.41e12,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE,
              RandDist="LogN")

### Joint Frailty model for clustered data ###

# here is a generated cluster (31 clusters of 13 subjects)
readmission <- transform(readmission,group=id\%\%31+1)

joi.clus <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(group)+
              num.id(id)+dukes+charlson+sex+chemo+terminal(death),
              formula.terminalEvent=~dukes+charlson+sex+chemo,
              data=readmission,Frailty=TRUE,joint=TRUE,recurrentAG=TRUE,
              n.knots=10,kappa1=2.11e+08,kappa2=9.53e+11)

}

}
\keyword{models}
