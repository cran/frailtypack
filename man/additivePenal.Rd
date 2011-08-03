\name{additivePenal}
\alias{additivePenal}
\Rdversion{1.1}
\title{Fit an Additive Frailty model using penalized likelihood estimation }

\description{

Fit an additive frailty model using penalized likelihood estimation.
The main issue in a meta-analysis study is how to take into account
the heterogeneity between trials and between the treatment effects across
 trials. Additive models are proportional hazards models with two correlated
 random trial effects that act either multiplicatively on the hazard function
 or in interaction with the treatment, which allows studying for instance meta-analysis or multicentric datasets.
 Right-censored data are allowed, but not the left-truncated data. A stratified analysis is possible (maximum number of strata: 2).
This approach is different from the shared gamma frailty models.  


In an additive model, the hazard function for the \eqn{j^{th}} subject in the \eqn{i^{th}} trial with random trial effect \eqn{u_i} as well as the random treatment-by-trial interaction \eqn{v_i} is:

\deqn{\lambda_{ij}(t|u_i,v_i)=\lambda_0(t)exp(u_i+v_iX_{ij1}+\sum_{k=1}^{p}\beta_kX_{ijk})}

\deqn{ u_i\sim\bold{\mathcal{N}}(0,\bold{\sigma^2}) \hspace{0.5cm} v_i\sim\bold{\mathcal{N}}(0,\bold{\tau^2})\hspace{0.5cm} \bold{cov}(u_i,v_i)=\bold{\rho\sigma\tau} }

where \eqn{\lambda_0(t)} is the baseline hazard function, \eqn{\beta_k} the fixed effect associated to the covariate \eqn{X_{ijk}} (k=1,..,p), \eqn{\beta_1} is the treatment effect and \eqn{X_{ij1}} the treatment variable. \eqn{\rho} is the corresponding correlation coefficient for the two frailty terms.

}


\usage{

additivePenal(formula, data, correlation = FALSE, recurrentAG =
                 FALSE, cross.validation = FALSE, n.knots, kappa1,
                 kappa2, maxit = 350)


}


\arguments{
  \item{formula}{ a formula object, with the response on the left of a \eqn{\texttildelow} operator, and the terms on the right. 
 The response must be a survival object as returned by the 'Surv' function like in survival
 package. The \code{slope()} function is required.}

  \item{data}{ a 'data.frame' in which to interpret the variables named in the
 'formula'.}

  \item{correlation}{ Logical value. Are the random effects correlated? If so,
 the correlation coefficient is estimated. The default is FALSE.}


  \item{recurrentAG}{Always FALSE for additive models (left-truncated data are not allowed).}

    \item{cross.validation}{Logical value. Is cross validation procedure used
 for estimating smoothing parameter?
                 If so a search of the smoothing parameter using cross
 validation is done, with kappa1 as the seed.
                 The cross validation is not implemented for two strata. The default is FALSE.}

  \item{n.knots}{ integer giving the number of knots to use. Value required.
  It corresponds to the (n.knots+2) splines functions for the approximation
 of the hazard or the survival functions. 
  Number of knots must be between 4 and 20.(See Note)}

  \item{kappa1}{ positive smoothing parameter. The coefficient kappa of the
 integral of the squared
          second derivative of hazard function in the fit (penalized log
 likelihood). To obtain an initial value for \code{kappa1} (or \code{kappa2}), a solution is to fit the corresponding shared frailty model using cross validation (See cross.validation).We advise the user 
          to identify several possible tuning parameters, note their
 defaults and look at the sensitivity 
          of the results to varying them. Value required.(See Note)}

 \item{kappa2}{Positive smoothing parameter for the second stratum when data are stratified. See kappa1.}
 
  \item{maxit}{maximum number of iterations for the Marquardt algorithm.
 Default is 350}

}


\details{
 
The estimated parameter are obtained by maximizing the penalized
 log-likelihood using the robust Marquardt algorithm (Marquardt,1963).
 The parameters are initialized with values obtained with Cox proportional
 hazard models. The iterations are stopped when the difference between two
 consecutive loglikelhoods was small \eqn{(<10^{-4})}, the estimated coefficients
 were stable (consecutive values \eqn{(<10^{-4})}, and the gradient small
 enough \eqn{(<10^{-6})}. To be sure of having a positive function at all stages of
 the algorithm, the spline coefficients were reparametrized to be positive
 at each stage. The variance space of the two random effects is reduced,
 so the variances are positive, and the correlation coefficient values are
 constrained to be  between -1 and 1. 
The marginal log-likelihood depends on integrations that are approximated
 by using the Laplace integration technique with a first order approximation.
 The smoothing parameter can be fixed or estimated by maximizing  likelihood
 cross-validation criterion. The usual squared Wald statistic was modified to
 a mixture of two \eqn{\chi^2} distribution to get significance test for the
 variance of the random effects.

\bold{INITIAL VALUES}

The splines and the regression coefficients are initialized to 0.1. An adjusted Cox model is fitted, it provides new initial values for the splines coefficients and the regression coefficients. The variances of the frailties are initialized to 0.1. Then an additive frailty model with independent frailties is fitted. At last, an additive frailty model with correlated frailties is fitted.


}

\value{
 An additive model or more generally an object of class 'additivePenal'.
 Methods defined for 'additivePenal' objects are provided for print, plot and
 summary. 


\item{b}{sequence of the corresponding estimation of the splines coefficients, the random effects variances and the regression coefficients.}
\item{call}{The code used for fitting the model}
\item{coef}{the coefficients of the linear predictor, which
          multiply the columns of the model matrix. }
\item{cov}{covariance between the two frailty terms \eqn{(\bold{cov}(u_i,v_i))}}
\item{cross.Val}{Logical value. Is cross validation procedure used for estimating the smoothing parameters?}
\item{correlation}{Logical value. Are the random effects correlated?}
\item{DoF}{degrees of freedom associated with the "kappa"}
\item{formula}{the formula part of the code used for the model}
\item{groups}{the maximum number of groups used in the fit}
\item{kappa}{ A vector with the smoothing parameters corresponding to each baseline function as components}
\item{lam}{matrix of hazard estimates at x1 times and confidence bands. }
\item{lam2}{the same value as lam for  the second stratum}
\item{loglikPenal}{the complete marginal penalized log-likelihood }

\item{n}{the number of observations used in the fit.}
\item{n.events}{the number of events observed in the fit}
\item{n.iter}{number of iterations needed to converge}
\item{n.knots }{number of knots for estimating the baseline function  }
\item{n.strat}{number of stratum}
\item{rho}{the corresponding correlation coefficient for the two frailty terms  }
\item{sigma2}{Variance for the random intercept (the random effect associated to the baseline hazard functon)  }
\item{surv}{matrix of baseline survival estimates at x1 times and confidence bands. }
\item{surv2}{the same value as surv for the  the second stratum}
\item{tau2}{Variance for the random slope (the random effect associated to the treatment effect across trials) }
%\item{type}{character string specifying the type of censoring. Possible values are %"right", "left", "counting", "interval", "interval2". The default is "right" or %"counting" depending on whether the 'time2' argument is absent (not interval-censored %data) or present (interval-censored data), respectively.}
 \item{varH}{the variance matrix of all parameters before positivity constraint transformation (Sigma2, Tau2, the regression coefficients and the spline coefficients). Thenafter, the delta method is needed to obtain the estimated variance parameters.}
 \item{varHIH}{the robust estimation of the variance matrix of all parameters (Sigma2, Tau2, the regression coefficients and the spline coefficients).}
 \item{varSigma2}{ The variance of the estimates of "sigma2"  }
 \item{varTau2}{ The variance of the estimates of "tau2"}
 \item{varcov}{ Variance of the estimates of "cov"  }
 \item{x1}{vector of times where both survival and hazard function for the recurrent events are estimated. By default seq(0,max(time),length=99), where time is the vector of survival times.}
 \item{x2}{see x1 value for the second stratum}

\item{npar}{number of parameters.}
\item{nvar}{number of vairables.}
\item{noVar}{indicator of explanatory variable.}
\item{LCV}{the approximated likelihood cross-validation criterion in the semi parametric case}
\item{n.knots.temp}{initial value of n.knots.}
}

\note{

"Kappa" (kappa1 and kappa2) and "n.knots" are the arguments that the user have to change if the fitted model does not converge. 
"n.knots" takes integer values between 4 and 20. But with n.knots=20, the model would take a long time to converge. So, usually, begin first with n.knots=7, and increase it step by step until it converges.
"Kappa" only takes positive values. So, choose a value for Kappa (for instance 10000), and if it does not converge, multiply or devide this value by 10 or 5 until it converges.

}



\references{ 

V. Rondeau, S. Michiels, B.Liquet, and J.P. Pignon (2008).
 Investigating trial and treatment heterogeneity in an individual patient data meta-analysis of survival data
 by mean of the maximum penalized likelihood approach.
\emph{Statistics in Medecine}, \bold{27}, 1894-1910.


}


\seealso{
 \code{\link{print.additivePenal}}, 
 \code{\link{plot.additivePenal}},
 \code{\link{summary.additivePenal}},
 \code{\link{cluster}},
\code{\link{slope}},
\code{\link{strata}},
\code{\link[survival]{Surv}}
 }


\examples{

 ### Additive model with 1 covariate ###

\dontrun{

  data(dataAdditive) 
  modAdd<-additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
                 correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=10000)


# It takes around 4 minutes to converge. Var1 is boolean as a treatment variable. #

  print(modAdd)
  summary(modAdd)
  plot(modAdd)

}


}

\keyword{file}
