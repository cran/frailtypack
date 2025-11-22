#' Fit a Weibull Illness-Death Model with Optional Shared Frailty
#'
#' @description Fit a three-state illness-death model (states: 0=Healthy, 1=Illness, 2=Death)
#' using Weibull baseline hazards for all transitions (0->1, 0->2, 1->2).
#' Allows for shared gamma frailty between the three transitions, acting
#' multiplicatively on the hazards within specified groups (clusters).
#' The model accommodates right-censored and left-truncated data.
#' The transition from the illness state to death (1->2) can be modeled using
#' either a Markov or a Semi-Markov assumption for the baseline hazards time scale.
#'
#' \if{html}{
#' \figure{IDSCHEME.png}{options: width="329"}
#' }
#' 
#' \if{latex}{
#' \figure{IDSCHEME.png}{options: width=8.7cm}
#' }
#'
#' @details Let \eqn{T_1} be the time to the non-terminal event (illness, 0->1) and
#' \eqn{T_2} be the time to the terminal event (death, 0->2 or 1->2).
#'
#' The transition intensities are defined as:
#' \deqn{
#'   \lambda_{01}(t) = \lim_{\Delta t \to 0^+} \frac{\mathbb{P}(t \leq T_1 \leq t + \Delta t \mid T_1 \geq t, T_2 \geq t)}{\Delta t}}
#'
#' \deqn{
#'   \lambda_{02}(t) = \lim_{\Delta t \to 0^+} \frac{\mathbb{P}(t \leq T_2 \leq t + \Delta t \mid T_1 \geq t, T_2 \geq t)}{\Delta t}}
#'
#' \deqn{
#'   \lambda_{12}(t \mid T_1 = s) = \lim_{\Delta t \to 0^+} \frac{\mathbb{P}(t \leq T_2 \leq t + \Delta t \mid T_1 = s, T_2 \geq t)}{\Delta t} \quad (0 < s < t)}
#'
#' A proportional hazards model with a shared frailty term \eqn{\omega_i} is assumed
#' for each transition within group \eqn{i}. For the \eqn{j^{th}} subject
#' (\eqn{j=1,...,n_i}) in the \eqn{i^{th}} group (\eqn{i=1,...,G}) the transition intensities are
#' defined as follows:
#' \deqn{
#'   \lambda_{01}^{ij}(t |\omega_i,X_{01}^{ij}) = \lambda_{0,01}(t) \omega_i \exp(\beta_1^{T} X_{01}^{ij})
#' }
#' \deqn{
#'   \lambda_{02}^{ij}(t |\omega_i,X_{02}^{ij}) = \lambda_{0,02}(t) \omega_i \exp(\beta_2^{T} X_{02}^{ij})
#' }
#' \deqn{
#'   \lambda_{12}^{ij}(t | T_1 = s, \omega_i, X_{12}^{ij}) = \lambda_{0,12}(t | T_1 = s) \omega_i \exp(\beta_3^{T} X_{12}^{ij}) \quad (0 < s < t)
#' }
#' \eqn{\omega_i} is the frailty term for the \eqn{i^{th}} group.

#' For subject-specific frailties, use \code{cluster(id)} where id is unique (\eqn{n_i=1}).

#' 

#' \eqn{\beta_1}, \eqn{\beta_2} and \eqn{\beta_3} are respectively the vectors 
#' of time fixed regression coefficients for the transitions 0->1, 0->2 and 1->2.

#' \eqn{X_{01}^{ij}}, \eqn{X_{02}^{ij}} and \eqn{X_{12}^{ij}} are respectively the vectors 
#' of time fixed covariates for the \eqn{j^{th}} subject in the \eqn{i^{th}} group for the 
#' transitions 0->1, 0->2 and 1->2.

#' \eqn{\lambda_{0,01}(.)}, \eqn{\lambda_{0,02}(.)} and \eqn{\lambda_{0,12}(.)} are respectively
#'  the baseline hazard functions for the transitions 0->1, 0->2 and 1->2.

#'


#'   The baseline hazard \eqn{\lambda_{0,12}(t | T_1 = s)} depends on the 'model' argument:
#'   \itemize{
#'   \item \bold{Markov model}: \eqn{\lambda_{0,12}(t | T_1 = s) = \lambda_{0,12}(t)}. 
#'   The risk depends on the time since origin.This model is suitable when the risk for
#'    the transition 1 -> 2 is primarily influenced by absolute time rather than the duration spent in state 1.
#'    
#'   \item \bold{Semi-Markov model}: \eqn{\lambda_{0,12}(t | T_1 = s) = \lambda_{0,12}(t - s)}. 
#'   The risk depends on the time since entering state 1 (sojourn time).
#'   This model is appropriate when the risk for the transition 1 -> 2  is
#'    more influenced by the time spent in state 1 rather than the time elapsed 
#'    since the initial starting point.
#'   }
#' The Weibull baseline hazard parameterization is:
#' \deqn{\lambda(t) = \frac{\gamma}{\lambda^\gamma} \cdot t^{\gamma - 1}} 


#' where \eqn{\lambda} is the scale parameter and \eqn{\gamma} the shape parameter
#'
#' 
#'
#' @usage frailtyIllnessDeath (formula, formula.terminalEvent, data, model = "Semi-Markov",
#' maxit = 300, init.B, init.Theta, init.hazard.weib,
#' LIMparam = 1e-3, LIMlogl = 1e-3, LIMderiv = 1e-3,
#' partialH, x01, x02, x12, print.info = FALSE,
#' print.result = TRUE, blinding = TRUE)
#' @param formula A formula object with the response on the left of a \eqn{\sim}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the 'Surv' function. Status should be 1 if event 0->1 occurred,
#' 0 otherwise. Covariates for transition 0->1 are specified here. Shared frailty
#' is specified using \code{cluster(group_variable)}. Left truncation can be
#' included via \code{Surv(time1, time2, status)}.
#' @param formula.terminalEvent A formula object. Response must be a \code{survival::Surv} object
#'  representing the time and status for the terminal event (death). Status should be 1 if
#'  transition to state 2 occurred (either from state 0 or state 1), 0 otherwise (censored). Covariates for
#'  transitions 0->2 and 1->2 are specified on the RHS (currently assumes the same covariates
#'  apply to both 0->2 and 1->2). 
#' 
#' \strong{Note on Left Truncation:} The illness-death model implemented assumes a
#' \strong{single entry time} per subject. This entry time, specified in \code{formula},
#' indicates that the subject must be in the initial state (state 0, having experienced
#' neither transition 0->1 nor 0->2) at that \code{entry_time} to be included
#' in the analysis. The entry time should \strong{not} be specified again here in
#' \code{formula.terminalEvent}.
#' @param data A 'data.frame' with the variables used in the formulas.
#' @param model Character string specifying the model for the 1->2 transition
#' baseline hazard. Allowed values: "Semi-Markov" (default) or "Markov".
#' @param maxit Maximum number of iterations for the Marquardt algorithm.
#' Default is 300.
#' @param init.B Optional. A vector of initial values for regression coefficients.
#'   Order is (\eqn{\beta_1},\eqn{\beta_2},\eqn{\beta_3}). 
#'   The total length must match the total number of regression coefficients. 
#'   If omitted, the regression coefficients are initialized from default values. 
#'   These defaults are obtained by fitting three independent Weibull proportional 
#'   hazards models (without frailty), one for each transition.
#' @param init.Theta Optional. Initial value for the frailty variance \eqn{\theta}.
#'   Default is 0.1. This parameter is only used if \code{cluster()} is present
#'   in \code{formula}.
#' @param init.hazard.weib Optional. A vector of initial values for the Weibull baseline
#'   hazard parameters. Must be of size 6, in the order: scale(0->1),
#'   shape(0->1), scale(0->2), shape(0->2), scale(1->2), shape(1->2).
#'   If omitted, the baseline hazard parameters are initialized from default values. 
#'   These defaults are obtained by fitting three independent Weibull proportional 
#'   hazards models (without frailty), one for each transition.
#' @param LIMparam Convergence threshold for the parameters based on the maximum
#'   absolute difference between successive iterations (\eqn{10^{-3}} by default).
#' @param LIMlogl Convergence threshold for the log-likelihood based on the absolute
#'   difference between successive iterations (\eqn{10^{-3}} by default).
#' @param LIMderiv Convergence threshold based on the relative distance to the optimum
#'   (related to gradient and Hessian) (\eqn{10^{-3}} by default). See Details.
#' @param x01 Optional. Numeric vector of time points at which to calculate baseline
#'   hazard and survival functions for transition 0->1. Defaults to a sequence of
#'   99 points from 0 to the maximum observed time for transition 0->1.
#' @param x02 Optional. Numeric vector of time points at which to calculate baseline
#'   hazard and survival functions for transition 0->2. Defaults to a sequence of
#'   99 points from 0 to the maximum observed time for transition 0->2.
#' @param x12 Optional. Numeric vector of time points at which to calculate baseline
#'   hazard and survival functions for transition 1->2. Defaults to a sequence of
#'   99 points from 0 to the maximum observed time for transition 1->2.
#' @param print.info Logical. If \code{TRUE}, prints information at each iteration
#'   of the optimization algorithm. Default is \code{FALSE}.
#' @param print.result Logical. If \code{TRUE}, prints a formatted summary of the
#'   results. Default is \code{TRUE}.
#' @param partialH Optional. Integer vector specifying the indices of parameters to
#'   exclude from the Hessian matrix when calculating the relative distance
#'   convergence criterion (\code{LIMderiv}). This is only considered if the first
#'   two criteria (\code{LIMparam}, \code{LIMlogl}) are met and the full Hessian
#'   is problematic (e.g., not invertible). Default is \code{NULL}.
#' @param blinding Logical. If \code{TRUE}, the algorithm attempts to continue even
#'   if the log-likelihood calculation produces non-finite values (e.g., \code{Inf},
#'   \code{NaN}) at some iteration. Setting to \code{FALSE} may cause the algorithm
#'   to stop earlier in such cases. Default is \code{TRUE}.
#'
#' @return An object of class 'frailtyIllnessDeath ' containing:
#' \describe{
#'   \item{b}{Vector of the estimated parameters. Order is:
#'   (scale(0->1), shape(0->1), scale(0->2), shape(0->2), scale(1->2), shape(1->2),
#'    \eqn{\hat {\theta}} (if frailty), \eqn{\hat{\beta}_1}, \eqn{\hat{\beta}_2}, \eqn{\hat{\beta}_3}).}
#'   \item{call}{The matched function call.}
#'   \item{coef}{Vector of estimated regression coefficients.}
#'   \item{loglik}{The marginal log-likelihood value at the final parameter estimates.}
#'   \item{grad}{Gradient vector of the log-likelihood at the final parameter estimates.}
#'   \item{n}{The number of subjects (observations) used in the fit.}
#'   \item{n.events}{Vector containing the number of observed events: count for 0->1, count for 0->2, count for 1->2 and count for censoring.}
#'   \item{n.iter}{Number of iterations.}
#'   \item{vcov}{Variance-covariance matrix for the parameters listed in \code{b}.}
#'   \item{npar}{Total number of estimated parameters.}
#'   \item{nvar}{Total number of regression coefficients.}
#'   \item{shape.weib}{Vector of estimated Weibull baseline shape parameters (shape(0->1),shape(0->2),shape(1->2)).}
#'   \item{scale.weib}{Vector of estimated Weibull baseline scale parameters (scale(0->1),scale(0->2),scale(1->2)).}
#'   \item{crit}{Convergence status code: 1=converged, 2=maximum iterations reached, 3=converged using partial Hessian, 4=the algorithm encountered a problem in the loglikelihood computation.}
#'   \item{Frailty}{Logical. \code{TRUE} if a model with shared frailty (\code{cluster(.)}) was fitted.}
#'   \item{beta_p.value}{Vector of p-values from Wald tests for the regression coefficients in \code{coef}.}
#'   \item{AIC}{Akaike Information Criterion, calculated as \eqn{AIC=\frac{1}{n}(np - l(.))}, where np is the number of parameters and l is the log-likelihood.}
#'   \item{x01}{Vector of time points used for calculating baseline functions for transition 0->1.}
#'   \item{x02}{Vector of time points used for calculating baseline functions for transition 0->2.}
#'   \item{x12}{Vector of time points (or sojourn times if Semi-Markov) used for calculating baseline functions for transition 1->2.}

#'   \item{lam01}{Matrix containing baseline hazard estimates and 95\% confidence intervals for transition 0->1  calculated at \code{x01}.}
#'   \item{lam02}{Matrix containing baseline hazard estimates and 95\% confidence intervals for transition 0->2 calculated at \code{x02}.}
#'   \item{lam12}{Matrix containing baseline hazard estimates and 95\% confidence intervals for transition 1->2 calculated at \code{x12}.}

#'   \item{surv01}{Matrix containing baseline survival estimates and 95\% confidence intervals for transition 0->1 calculated at \code{x01}.}
#'   \item{surv02}{Matrix containing baseline survival estimates and 95\% confidence intervals for transition 0->2 calculated at \code{x01}.}
#'   \item{surv12}{Matrix containing baseline survival estimates and 95\% confidence intervals for transition 0->2 calculated at \code{x12}.}

#'   \item{median.01}{Matrix containing the estimated median baseline survival time and its 95\% confidence interval for transition 0->1.}
#'   \item{median.02}{Matrix containing the estimated median baseline survival time and its 95\% confidence interval for transition 0->2.}
#'   \item{median.12}{Matrix containing the estimated median baseline survival time and its 95\% confidence interval for transition 1->2.}

#'   \item{linear.pred01}{Vector of linear predictors calculated for transition 0->1. For non-frailty models, this is \eqn{\hat{\beta}_1^{t}X_{01}}. For frailty models, it includes the estimated log-frailty: \eqn{\hat{\beta}_1^{t}X_{01} + \log(\hat{\omega}_i).}}
#'   \item{linear.pred02}{Vector of linear predictors calculated for transition 0->2. For non-frailty models, this is \eqn{\hat{\beta}_2^{t}X_{02}}. For frailty models, it includes the estimated log-frailty: \eqn{\hat{\beta}_2^{t}X_{02} + \log(\hat{\omega}_i).}}
#'   \item{linear.pred12}{Vector of linear predictors calculated for transition 1->2. For non-frailty models, this is \eqn{\hat{\beta}_3^{t}X_{12}}. For frailty models, it includes the estimated log-frailty: \eqn{\hat{\beta}_3^{t}X_{12} + \log(\hat{\omega}_i).}}

#'   \item{names.factor.01}{Character vector identifying factor covariates included for transition 0->1.}
#'   \item{names.factor.02}{Character vector identifying factor covariates included for transition 0->2.}
#'   \item{names.factor.12}{Character vector identifying factor covariates included for transition 1->2.}

#'   \item{global_chisq.01}{Vector containing the chi-squared statistics for global Wald tests of factor variables for transition 0->1.}
#'   \item{global_chisq.02}{Vector containing the chi-squared statistics for global Wald tests of factor variables for transition 0->2.}
#'   \item{global_chisq.12}{Vector containing the chi-squared statistics for global Wald tests of factor variables for transition 1->2.}

#'   \item{dof_chisq.01}{Vector containing the degrees of freedom for the global Wald tests for transition 0->1.}
#'   \item{dof_chisq.02}{Vector containing the degrees of freedom for the global Wald tests for transition 0->2.}
#'   \item{dof_chisq.12}{Vector containing the degrees of freedom for the global Wald tests for transition 1->2.}

#'   \item{p.global_chisq.01}{Vector containing the p-values for the global Wald tests for transition 0->1.}
#'   \item{p.global_chisq.02}{Vector containing the p-values for the global Wald tests for transition 0->2.}
#'   \item{p.global_chisq.12}{Vector containing the p-values for the global Wald tests for transition 1->2.}

#'   \item{global_chisq.test.01}{Indicator (0/1) whether any global factor tests were performed for transition 0->1.}
#'   \item{global_chisq.test.02}{Indicator (0/1) whether any global factor tests were performed for transition 0->2.}
#'   \item{global_chisq.test.12}{Indicator (0/1) whether any global factor tests were performed for transition 1->2.}

#' }
#' If \code{Frailty} is \code{TRUE}, the following components related to frailty are also included:
#' \describe{
#'   \item{groups}{The number of unique groups specified by \code{cluster(.)}.}
#'   \item{theta}{The estimated variance (\eqn{\hat{\theta}}) of the Gamma frailty distribution.}
#'   \item{theta_p.value}{The p-value from a Wald test for the null hypothesis \eqn{H_0: \theta=0}.}
#'   \item{VarTheta}{The estimated variance of the frailty variance estimator: \eqn{\hat{Var}(\hat{\theta})}.}
#'   \item{frailty.pred}{Vector containing the empirical Bayes predictions of the frailty term for each group.}
#'   \item{frailty.var}{Vector containing the variances of the empirical Bayes frailty predictions.}
#'   \item{frailty.sd}{Vector containing the standard errors of the empirical Bayes frailty predictions.}
#' }
#' @note The optimization uses the robust Marquardt algorithm (Marquardt, 1963),
#' combining Newton-Raphson and steepest descent steps. Iterations stop when
#' criteria \code{LIMparam}, \code{LIMlogl}, and \code{LIMderiv} are all met.
#' Confidence bands for the baseline hazard and baseline survival functions were 
#' computed using a Monte Carlo simulation approach based on the estimated Weibull 
#' parameters. A sample of size 1000 was drawn from the joint distribution of the 
#' shape and scale parameters. For each sampled parameter set, the baseline functions 
#' were evaluated over the time grid defined by the vector \code{x0}. Pointwise 95% 
#' confidence bands were then obtained by computing the 2.5th and 97.5th percentiles 
#' of the simulated values at each time point for each baseline function.

#'
#' @references
#' Lee, C., Gilsanz, P., & Haneuse, S. (2021). Fitting a shared frailty
#' illness-death model to left-truncated semi-competing risks data
#' to examine the impact of education level on incident dementia.
#' \emph{BMC Medical Research Methodology}, 21(1), 1-13.
#'
#' Marquardt, D. W. (1963). An algorithm for least-squares estimation of nonlinear
#' parameters. \emph{SIAM Journal on Applied Mathematics}, 11(2), 431-441.
#'
#' Liquet, B., Timsit, J. F., & Rondeau, V. (2012). Investigating hospital
#' heterogeneity with a multi-state frailty model: application to nosocomial
#' pneumonia disease in intensive care units. \emph{BMC Medical Research
#' Methodology}, 12(1), 1-14.

#' @importFrom dplyr %>% filter  
#' @examples

#' \donttest{

#'  
#' data(readmission2)
#'
#' ##-- Fitting models --##
#'
#' #-- Semi-Markovian Weibull Illness-Death model with 
#' #   group frailty shared between transitions --#
#'
#' ModIllnessDeath_SemiMarkov_Group <- frailtyIllnessDeath(
#'   formula = Surv(observed_disease_time, disease_status) ~ cluster(group) + sex,
#'   formula.terminalEvent = Surv(observed_death_time, death_status) ~ sex,
#'   data        = readmission2,
#'   model       = "Semi-Markov",
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' #-- Markovian Weibull Illness-Death model with subject-specific 
#' #   frailty shared between transitions --#
#'
#' ModIllnessDeath_Markov_Subject <- frailtyIllnessDeath(
#'   formula = Surv(observed_disease_time, disease_status) ~ cluster(id) + sex,
#'   formula.terminalEvent = Surv(observed_death_time, death_status) ~ sex,
#'   data        = readmission2,
#'   model       = "Markov",
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' #--- Semi-Markovian Weibull Illness-Death model with a factor and 
#' #    no covariates for the non-terminal event ---#
#'
#' ModIllnessDeath_SemiMarkov_NoCov_Factor <- frailtyIllnessDeath(
#'   formula = Surv(observed_disease_time, disease_status) ~ 1,
#'   formula.terminalEvent = Surv(observed_death_time, death_status) ~ factor(dukes),
#'   data        = readmission2,
#'   model       = "Semi-Markov",
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' #--- Semi-Markovian Weibull Illness-Death model with left truncation ---#
#'
#' data(Paq810)
#'
#' ModIllnessDeath_SemiMarkov_LeftTrunc <- frailtyIllnessDeath(
#'   formula = Surv(e, r, dementia) ~ gender + certif,
#'   formula.terminalEvent = Surv(t, death) ~ certif,
#'   data        = Paq810,
#'   model       = "Semi-Markov",
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' #--- Markovian Weibull Illness-Death model with left truncation ---#
#'
#' ModIllnessDeath_Markov_LeftTrunc <- frailtyIllnessDeath(
#'   formula = Surv(e, r, dementia) ~ gender + certif,
#'   formula.terminalEvent = Surv(t, death) ~ certif,
#'   data        = Paq810,
#'   model       = "Markov",
#'   print.info  = FALSE,
#'   maxit       = 100
#' ) }
#'
#' @importFrom dplyr %>% filter      
#' @importFrom tidyr drop_na




















#' @export
"frailtyIllnessDeath"  <- function(formula, 
                          formula.terminalEvent, data, model="Semi-Markov", maxit = 300,
                          init.B, init.Theta,init.hazard.weib,
                          LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, partialH,
                          x01, x02, x12, print.info=FALSE, print.result=TRUE, blinding = TRUE
)
{	
  
  
  ################
  ## Log-likelihood illness-death, shared frailty, left-truncated data
  ## Weibull baseline hazards
  logLike.weibull.SCR.SM.LT.SANS.FRAILTY.SEMI.MARKOV <- function(para, y1, y2, delta1, delta2,l=l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL)
  {
    ##
    kappa1    <- exp(para[1])
    alpha1 <- exp(para[2])
    kappa2    <- exp(para[3])
    alpha2 <- exp(para[4])
    kappa3    <- exp(para[5])
    alpha3 <- exp(para[6])
    
    nP.0 <- 6
    nP.1 <- ncol(Xmat1)
    nP.2 <- ncol(Xmat2)
    nP.3 <- ncol(Xmat3)
    ##
    eta.1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
    eta.2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
    eta.3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
    ##
    
    ##
    
    
    type1 <- as.numeric(delta1 == 1 & delta2 == 1 )
    type2 <- as.numeric(delta1 == 1 & delta2 == 0 )
    type3 <- as.numeric(delta1 == 0 & delta2 == 1 )
    type4 <- as.numeric(delta1 == 0 & delta2 == 0 )
    
    
    log.h1star.y1 <- log(alpha1) -alpha1*log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
    log.h2star.y1 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
    log.h2star.y2 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
    log.h3star.y2 <- log(alpha3) -alpha3* log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
    ##
    q.y1 <- (y1/kappa1)^alpha1 * exp(eta.1) + (y1/kappa2)^alpha2 * exp(eta.2)
    q.l <- (l/kappa1)^alpha1 * exp(eta.1) + (l/kappa2)^alpha2 * exp(eta.2)
    ##
    w.y1.y2 <- ((y2-y1)/kappa3)^alpha3 * exp(eta.3)
    ##
    k1 <- w.y1.y2
    k2.y1 <- q.y1 
    ##
    
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1) + q.l
    logLike2 <- log.h1star.y1 - (k1 + k2.y1) +  q.l
    logLike3 <- log.h2star.y1 - k2.y1 +  q.l
    logLike4 <- - k2.y1  +  q.l
    
    
    
    loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) 
    ##
    return(loglh)
    
    
  }
  
  
  
  ################
  ## Log-likelihood illness-death, shared frailty, left-truncated data
  ## Weibull baseline hazards
  logLike.group.weibull.SCR.SM.LT.FRAILTY.SEMI.MARKOV <- function(para, y1, y2, delta1, delta2, l,group,data, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL)
  {
    ##
    kappa1    <- exp(para[1])
    alpha1 <- exp(para[2])
    kappa2    <- exp(para[3])
    alpha2 <- exp(para[4])
    kappa3    <- exp(para[5])
    alpha3 <- exp(para[6])
    
    theta    <- ((para[7])^2)
    thetaInv <- 1 / theta
    
    ##
    nP.0 <- 7
    nP.1 <- ncol(Xmat1)
    nP.2 <- ncol(Xmat2)
    nP.3 <- ncol(Xmat3)
    ##
    eta.1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
    eta.2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
    eta.3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
    ##
    
    ##
    
    
    log.h1star.y1 <- log(alpha1) -alpha1*log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
    log.h2star.y1 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
    log.h3star.y2 <- log(alpha3) -alpha3* log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
    
    
    q.y1 <- (y1/kappa1)^alpha1 * exp(eta.1) + (y1/kappa2)^alpha2 * exp(eta.2)
    k1 <- delta1*((y2-y1)/kappa3)^alpha3 * exp(eta.3)
    
    q.l <- (l/kappa1)^alpha1 * exp(eta.1) + (l/kappa2)^alpha2 * exp(eta.2)
    
    
    loglike <- 0
    check <- c()
    for(i in unique(group)){
      
      ni <- which(data$group == i)
      mi1 <- length(which(delta1[ni]==1)) ## sum(delta1[ni])
      mi2 <- length(which(delta2[ni]==1) )
      
      for(j in ni){
        loglike <- loglike + delta1[j]*log.h1star.y1[j]+ delta2[j]*(1-delta1[j])*log.h2star.y1[j]
        if(delta1[j]==1 && delta2[j]==1){
          loglike <- loglike + (delta2[j]*delta1[j])*log.h3star.y2[j]
        }
      }
      
      loglike <- loglike - thetaInv*log(theta)
      
      if(mi1+mi2>0){
        for(k in 1:(mi1+mi2)){
          loglike <- loglike + log(mi1+mi2+thetaInv-k)
          
        }
      }
      
      
      v=0
      for(j in ni){
        v <- v + (q.y1[j]) + (k1[j])
      }
      
      loglike <- loglike - (mi1+mi2+thetaInv)*log(v+thetaInv)
      
      
      ### LEFT TRUNCATION TERM
      
      lt_sum = 0
      
      for (j in ni){
        lt_sum <- lt_sum + (q.l[j])
      }
      
      loglike <- loglike+ thetaInv*log(1+theta*lt_sum)
      
    }
    
    
    
    return(loglike)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###MARKOV 
  
  
  
  ################
  ## Log-likelihood illness-death, shared frailty, left-truncated data
  ## Weibull baseline hazards
  logLike.weibull.SCR.SM.LT.SANS.FRAILTY.MARKOV <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL)
  {
    ##
    kappa1    <- exp(para[1])
    alpha1 <- exp(para[2])
    kappa2    <- exp(para[3])
    alpha2 <- exp(para[4])
    kappa3    <- exp(para[5])
    alpha3 <- exp(para[6])
    
    ##
    nP.0 <- 6
    nP.1 <- ncol(Xmat1)
    nP.2 <- ncol(Xmat2)
    nP.3 <- ncol(Xmat3)
    ##
    eta.1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
    eta.2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
    eta.3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
    ##
    
    ##
    
    
    type1 <- as.numeric(delta1 == 1 & delta2 == 1 )
    type2 <- as.numeric(delta1 == 1 & delta2 == 0 )
    type3 <- as.numeric(delta1 == 0 & delta2 == 1 )
    type4 <- as.numeric(delta1 == 0 & delta2 == 0 )
    
    
    log.h1star.y1 <- log(alpha1) -alpha1*log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
    log.h2star.y1 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
    log.h2star.y2 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
    log.h3star.y2 <- log(alpha3) -alpha3* log(kappa3) + (alpha3 - 1) * log(y2) + eta.3
    ##
    q.y1 <- (y1/kappa1)^alpha1 * exp(eta.1) + (y1/kappa2)^alpha2 * exp(eta.2)
    q.y2 <- (y2/kappa1)^alpha1 * exp(eta.1) + (y2/kappa2)^alpha2 * exp(eta.2)
    q.l <- (l/kappa1)^alpha1 * exp(eta.1) + (l/kappa2)^alpha2 * exp(eta.2)
    ##
    w.y1.y2 <- (((y2^alpha3)-(y1^alpha3))/(kappa3^alpha3)) * exp(eta.3)
    ##
    k1 <- w.y1.y2
    k2.y1 <- q.y1 
    k2.y2 <- q.y2 
    ##
    
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1) +  q.l
    logLike2 <- log.h1star.y1 - (k1 + k2.y1)  +  q.l
    logLike3 <- log.h2star.y1 - k2.y1 +  q.l
    logLike4 <- - k2.y1  +  q.l
    
    
    
    loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) 
    ##
    
    
    return(loglh)
  }
  
  
  
  ################
  ## Log-likelihood illness-death, shared frailty, left-truncated data
  ## Weibull baseline hazards
  logLike.group.weibull.SCR.SM.LT.FRAILTY.MARKOV <- function(para, y1, y2, delta1, delta2, l,group,data, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL)
  {
    ##
    kappa1    <- exp(para[1])
    alpha1 <- exp(para[2])
    kappa2    <- exp(para[3])
    alpha2 <- exp(para[4])
    kappa3    <- exp(para[5])
    alpha3 <- exp(para[6])
    
    theta    <- ((para[7])^2)
    thetaInv <- 1 / theta
    
    ##
    nP.0 <- 7
    nP.1 <- ncol(Xmat1)
    nP.2 <- ncol(Xmat2)
    nP.3 <- ncol(Xmat3)
    ##
    eta.1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
    eta.2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
    eta.3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
    ##
    
    ##
    
    
    log.h1star.y1 <- log(alpha1) -alpha1*log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
    log.h2star.y1 <- log(alpha2) -alpha2* log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
    log.h3star.y2 <- log(alpha3) -alpha3* log(kappa3) + (alpha3 - 1) * log(y2) + eta.3
    
    
    q.y1 <- (y1/kappa1)^alpha1 * exp(eta.1) + (y1/kappa2)^alpha2 * exp(eta.2)
    k1 <- delta1*(((y2^alpha3)-(y1^alpha3))/(kappa3^alpha3)) * exp(eta.3)
    
    q.l <- (l/kappa1)^alpha1 * exp(eta.1) + (l/kappa2)^alpha2 * exp(eta.2)
    
    
    loglike <- 0
    check <- c()
    for(i in unique(group)){
      
      ni <- which(data$group == i)
      mi1 <- length(which(delta1[ni]==1)) ## sum(delta1[ni])
      mi2 <- length(which(delta2[ni]==1) )
      
      for(j in ni){
        loglike <- loglike + delta1[j]*log.h1star.y1[j]+ delta2[j]*(1-delta1[j])*log.h2star.y1[j]
        if(delta1[j]==1 && delta2[j]==1){
          loglike <- loglike + (delta2[j]*delta1[j])*log.h3star.y2[j]
        }
      }
      
      loglike <- loglike - thetaInv*log(theta)
      
      if(mi1+mi2>0){
        for(k in 1:(mi1+mi2)){
          loglike <- loglike + log(mi1+mi2+thetaInv-k)
          
        }
      }
      
      
      v=0
      for(j in ni){
        v <- v + (q.y1[j]) + (k1[j])
      }
      
      loglike <- loglike - (mi1+mi2+thetaInv)*log(v+thetaInv)
      
      
      ### LEFT TRUNCATION TERM
      
      lt_sum = 0
      
      for (j in ni){
        lt_sum <- lt_sum + (q.l[j])
      }
      
      loglike <- loglike+ thetaInv*log(1+theta*lt_sum)
      
    }
    
    
    
    
    return(loglike)
  }
  
  
  
  
  
  
  # Function to calculate confidence bands for a set of times
  
  hazard_and_confidence_bands_illdeath <- function(t_vec,scale,shape,varcov) {

    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    
    ##DELTA METHOD TO OBTAIN THE VARCOV MATRIX FOR LOG OF SHAPE AND SCALE
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)
    
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline hazard estimate","Lower", "Upper")    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      haz_pert <- numeric(2000)
      for (k in 1:2000) {
        haz_pert[k] <- (Sample_scale_shape[k,2]/ (Sample_scale_shape[k,1]^Sample_scale_shape[k,2])) * (t^(Sample_scale_shape[k,2] - 1))
      }
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      result[t_idx, 3] <- quantile(haz_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(haz_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  survival_and_confidence_bands_illdeath <- function(t_vec,scale,shape,varcov) {
    
    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)  
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline survival estimate","Lower", "Upper")
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      surv_pert <- numeric(2000)
      for (k in 1:2000) {
        surv_pert[k] <- exp(-(t/Sample_scale_shape[k,1])^Sample_scale_shape[k, 2])  }
      
      result[t_idx,1] <- t
      result[t_idx,2] <- exp(-(t/scale)^shape)
      
      result[t_idx, 3] <- quantile(surv_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(surv_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  
  
  

  
  
  
  base_hazard <- function(t_vec,scale,shape) {
    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline hazard estimate")
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      
    }
    
    return(result)
  }
  
  
  
  
  
  su <- function(t_vec,scale,shape) {

    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline survival estimate")
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      result[t_idx,1] <- t
      result[t_idx,2] <- exp(-(t/scale)^shape)
      
      
    }
    
    return(result)
  }
  
  
  
  
  
  
  
  minmin <- function(y, x) {
    tolerance <- .Machine$double.eps^.5  
    keep <- (!is.na(y) & y <(.5 + tolerance))
    if (!any(keep)) NA
    else {
      x <- x[keep]
      y <- y[keep]
      if (abs(y[1]-.5) <tolerance  && any(y< y[1])) 
        (x[1] + x[min(which(y<y[1]))])/2
      else x[1]
    }
  }
  
  
  
  
  
  direct <- c()
  
  if(!missing(print.result)){
    if (!is.logical(print.result) || length(print.result) != 1) {
      stop("'print.result' should be a logical value (TRUE or FALSE).")
    }
  }
  
  

  
  if(print.result==TRUE){
    direct=TRUE
  }
  
  if (print.result==FALSE){
    direct=FALSE
  }

  if(!missing(print.info)){
    if (!is.logical(print.info) || length(print.info) != 1) {
      stop("'print.info' should be a logical value (TRUE or FALSE).")
    }
  }
  

  
  
  
  if(!missing(blinding)){
    if (!is.logical(blinding) || length(blinding) != 1) {
      stop("'blinding' should be a logical value (TRUE or FALSE).")
    }
  }
  
  
  
  
  
  if(!missing(LIMparam)){
    if (!is.numeric(LIMparam) || length(LIMparam) != 1 || LIMparam <= 0) {
      stop("'LIMparam' should be a positive numeric value.")
    }
  }
  
  if(!missing(LIMderiv)){
    if (!is.numeric(LIMderiv) || length(LIMderiv) != 1 || LIMderiv <= 0) {
      stop("'LIMderiv' should be a positive numeric value.")
    }
  }
  
  if(!missing(LIMlogl)){
    if (!is.numeric(LIMlogl) || length(LIMlogl) != 1 || LIMlogl <= 0) {
      stop("'LIMlogl' should be a positive numeric value.")
    }
  }
  
  if (!missing(maxit) && (!is.numeric(maxit) || length(maxit) != 1 || maxit <= 0 || maxit != round(maxit))) {
    stop("'maxit' should be a positive integer.")
  }
  
  if (!is.data.frame(data)) {
    stop("'data' should be a data frame.")
  }
  
  
  
  
  
  valid_models <- c("Markov", "Semi-Markov")
  if(!missing(model)){
    if (!(model %in% valid_models)) {
      stop("'model' should be one of: 'Markov', or 'Semi-Markov'.")
    }
    
  }
  
  
  
  
  if(!missing(init.hazard.weib)){
    
    list_init.hazard.weib <- list()
    for(i in 1:length(init.hazard.weib)) {
      if(is.na(suppressWarnings({as.numeric(init.hazard.weib[i])}))){
        list_init.hazard.weib <- append(list_init.hazard.weib,as.character(init.hazard.weib[i]))
      }
      else{
        list_init.hazard.weib <- append(list_init.hazard.weib,as.numeric(init.hazard.weib[i]))
      }
    }
    
    init.hazard.weib <- list_init.hazard.weib
    
    
    if(any(sapply(init.hazard.weib, function(x) !(is.numeric(x) && x>0) ))) {
      stop("init.hazard.weib should be a vector of positive real numbers.")
    }
    
    if(length(init.hazard.weib)!=6){
      stop("Wrong number of baseline hazards coefficients in init.hazard.weib, this vector should be of 
           size 6 (a shape and a scale for each of the 3 transitions)")
    }
  } 
  
  
  
  
  
  
  if(!missing(init.Theta)){
    if (!is.numeric(init.Theta) || length(init.Theta) != 1 || init.Theta < 0) {
      stop("'init.Theta' should be a positive numeric value.")
    }
  }
  
  

  if (!inherits(formula, "formula") || !inherits(formula.terminalEvent, "formula")) {
    stop("'formula' and 'formula.terminalEvent' should be formulas.")
  }
  
  
  terms_object1 <- terms(formula)
  term_labels1 <- attr(terms_object1, "term.labels")
  term_labels1 <- term_labels1[!grepl("^cluster\\(.*\\)", term_labels1)]  # remove cluster()
  
  covariates_for_check <- gsub("^factor\\((.*)\\)$", "\\1", term_labels1)
  covariates_for_check <- unlist(lapply(covariates_for_check, function(term) {
    if (grepl("^I\\(.*\\)$", term)) {
      inner <- gsub("^I\\((.*)\\)$", "\\1", term)
      all.vars(parse(text = inner))
    } else {
      term
    }
  }))
  
  missing_covariates1 <- setdiff(as.character(covariates_for_check), names(data))
  if (length(missing_covariates1) > 0) {
    stop(paste(
      "Some covariates in formula were not found in the data:",
      paste(missing_covariates1, collapse = ", ")
    ))
  }
  
  
  terms_object2 <- terms(formula.terminalEvent)
  term_labels2 <- attr(terms_object2, "term.labels")
  term_labels2 <- term_labels2[!grepl("^cluster\\(.*\\)", term_labels2)]  # remove cluster()
  
  covariates_for_check2 <- gsub("^factor\\((.*)\\)$", "\\1", term_labels2)
  covariates_for_check2 <- unlist(lapply(covariates_for_check2, function(term) {
    if (grepl("^I\\(.*\\)$", term)) {
      inner <- gsub("^I\\((.*)\\)$", "\\1", term)
      all.vars(parse(text = inner))
    } else {
      term
    }
  }))
  
  missing_covariates2 <- setdiff(as.character(covariates_for_check2), names(data))
  if (length(missing_covariates2) > 0) {
    stop(paste(
      "Some covariates in formula.terminalEvent were not found in the data:",
      paste(missing_covariates2, collapse = ", ")
    ))
  }
  
  if (any(grepl("^cluster\\(.*\\)", attr(terms(formula.terminalEvent), "term.labels")))) {
    stop("Error in formula.terminalEvent: For a SHARED frailty model, the cluster term must be specified only in 'formula'.")
  }
  
  
  
  
  
  
  if(!missing(init.B)){
    
    list_init.B <- list()
    for(i in 1:length(init.B)) {
      if(is.na(suppressWarnings({as.numeric(init.B[i])}))){
        list_init.B <- append(list_init.B,as.character(init.B[i]))
      }
      else{
        list_init.B <- append(list_init.B,as.numeric(init.B[i]))
      }
    }
    
    init.B <- list_init.B
    
    
    if(any(sapply(init.B, function(x) !is.numeric(x) ))) {
      stop("init.B should be a vector of real numbers.")
    }
    
    if(length(init.B)!=(length(covariates_for_check)+2*length(covariates_for_check2))){
      stop("Wrong number of regression coefficients in init.B")
    }
  }
  
  
  
  
  
  
  
  Frailty=0
  
  if(any(grepl("cluster\\(.*\\)", deparse(formula)))){
    Frailty=1
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  
  rhs_expr <- formula[[3]]
  
  rhs_terms <- list()
  
  if (is.call(rhs_expr)) {
    rhs_terms <- unlist(strsplit(deparse(rhs_expr), split = " +"))  # Split on "+" symbol
  } else {
    rhs_terms <- list(deparse(rhs_expr))
  }
  
  cluster_term <- NULL
  if (is.call(rhs_expr) && ("+" %in% rhs_terms)) {
    rhs_terms <- as.list(rhs_terms)
    for (term in rhs_terms) {
      if ( grepl("^\\w+\\(.*\\)$", term) && sub("\\((.*?)\\)", "", term) == "cluster") {
        cluster_term <- term
        break
      }
    }
  }
  
  
  if(!("+" %in% rhs_terms)){
    if ( grepl("^\\w+\\(.*\\)$", rhs_terms) &&  sub("\\((.*?)\\)", "", rhs_terms) == "cluster") {
      cluster_term <- rhs_terms
    }
  }
  if (!is.null(cluster_term)) {
    cluster_variable <- sub("cluster\\((.*?)\\)", "\\1", deparse(cluster_term))
    group <- eval(parse(text = paste0("data$", cluster_variable)))
    data$group <- group
  }
  
  if(!missing(init.Theta)){
    
    if(is.null(cluster_term)){ 
      stop("init.Theta does not exist in a Weibull Illness Death model")
    }
    
  }
  
  
  
  
  
  
  
  lhs <- formula[[2]]
  
  
  if (inherits(lhs, "call") && lhs[[1]] == as.symbol("Surv")) {
    surv_terms1 <- as.character(lhs[-1])  
    
    
    
    surv_trunc1 <- surv_terms1
    
    if(length(surv_terms1)==2){
      if(!grepl("\\$",surv_terms1[1])){
        surv_terms1[1] <- paste0("data$",surv_terms1[1])
      }
      
      if(!grepl("\\$",surv_terms1[2])){
        surv_terms1[2] <- paste0("data$",surv_terms1[2])
      }
    }
    
    
    
    
    if(length(surv_terms1)==3){
      if(!grepl("\\$",surv_terms1[1])){
        surv_terms1[1] <- paste0("data$",surv_terms1[1])
      }
      
      if(!grepl("\\$",surv_terms1[2])){
        surv_terms1[2] <- paste0("data$",surv_terms1[2])
      }
      
      if(!grepl("\\$",surv_terms1[3])){
        surv_terms1[3] <- paste0("data$",surv_terms1[3])
      }
      
      
    }
    
    
    formula11 <- paste("Surv(", paste(surv_terms1, collapse = ", "), ")")
    formula11 <- eval(parse(text = formula11))
    
  }
  
  
  
  if(!inherits(lhs, "call")) {
    if(inherits(eval(lhs), "Surv")){
      formula11 <- eval(lhs)  
      surv_trunc1 <- attr(eval(lhs), "dimnames")[[2]]
      surv_terms1 <- surv_trunc1
    }
  }
  
  if (!inherits(formula11, "Surv")) {
    stop("The left hand side of 'formula' must be an object of type 'Surv'.")
  }
  
  
  lhs <- formula.terminalEvent[[2]]
  
  
  if (inherits(lhs, "call") && lhs[[1]] == as.symbol("Surv")) {
    surv_terms2 <- as.character(lhs[-1])  
    
    
    if(length(surv_terms2)==3){
      
      stop("Left truncation time is only specified in 'formula'.")
    }
    
    
    if(length(surv_terms2)==2){
      if(!grepl("\\$",surv_terms2[1])){
        surv_terms2[1] <- paste0("data$",surv_terms2[1])
      }
      
      if(!grepl("\\$",surv_terms2[2])){
        surv_terms2[2] <- paste0("data$",surv_terms2[2])
      }
    }
    
    
    
    
  
    
    
    formula22 <- paste("Surv(", paste(surv_terms2, collapse = ", "), ")")
    formula22 <- eval(parse(text = formula22))
  }
  
  
  
  if(!inherits(lhs, "call")) {
    if(inherits(eval(lhs), "Surv") ){
      formula22 <- eval(lhs)  
      surv_terms2 <- attr(eval(lhs), "dimnames")[[2]]
      
    }
    
    
    if(length(surv_terms2)==3){
      
      stop("Left truncation time is only specified in 'formula'.")
    }
  }
  
  
  if (!inherits(formula22, "Surv")) {
    stop("The left hand side of 'formula.terminalEvent' must be an object of type 'Surv'.")
  }
  
  
  

    vars_to_keep <- unique(c(
    covariates_for_check,      
    covariates_for_check2,     
    surv_terms1,               
    surv_terms2               
  ))
  
  if (!is.null(cluster_term)) {
    cluster_variable <- sub("cluster\\((.*?)\\)", "\\1", cluster_term)
    vars_to_keep <- unique(c(vars_to_keep, cluster_variable))
  }
  
  vars_to_keep <- vars_to_keep[vars_to_keep %in% names(data)]
  
  data <- data[complete.cases(data[, vars_to_keep, drop = FALSE]), , drop = FALSE]
  
  
  
  
  
  
  
  if(length(surv_trunc1)==3){
    y1     <- as.vector(formula11[,"stop"])
    delta1 <- as.vector(formula11[,"status"])
    l <-  as.vector(formula11[,"start"])
  }
  
  if(length(surv_trunc1)==2){
    y1     <- as.vector(formula11[,"time"])
    delta1 <- as.vector(formula11[,"status"])
    l <-  as.vector(rep(0,nrow(data)))
  }
  
  
  
  
  
  y2     <- as.vector(formula22[,"time"])
  delta2 <- as.vector(formula22[,"status"])
  
  
  if(missing(x01)){
    x01=round(seq(0,max(y1),length=99),3)
  }
  
  if(missing(x02)){
    x02=round(seq(0,max(y2[which(delta1==0)]),length=99),3)
  }
  
  if(missing(x12)){
    x12=round(seq(0,max(y2[which(delta1==1)]-y1[which(delta1==1)]),length=99),3)
  }
  
  
  
  if(!missing(x01)){
    if(!all(sapply(x01, function(x) is.numeric(x) && x >= 0))){
      stop("The vector of times for the transition 0->1 'x01' should be a vector of positive real numbers.")
    }
  }
  
  
  if(!missing(x02)){
    if(!all(sapply(x02, function(x) is.numeric(x) && x >= 0))){
      stop("The vector of times for the transition 0->2 'x02' should be a vector of positive real numbers.")
    }
  }
  
  if(!missing(x12)){
    if(!all(sapply(x12, function(x) is.numeric(x) && x >= 0))){
      stop("The vector of times for the transition 1->2 'x12' should be a vector of positive real numbers.")
    }
  }
  
  data_init <- data.frame(y1=y1,delta1=delta1,
                          y2=y2,delta2=delta2,
                          l=l)
  

  ###XMAT1
  
  covfactors1 <- c()
  terms_object1 <- terms(formula)
  covariates1 <- attr(terms_object1, "term.labels")
  covariates_without_cluster1 <- covariates1[!grepl("^cluster\\(.*\\)", covariates1)]
  
  
  
  if(length(covariates_without_cluster1)>0){
    
    Xmat1 <- NULL
    
    
    for(i in 1:length(covariates_without_cluster1)){
      
      columnname <- c()
      
      
      if(grepl("\\$",covariates_without_cluster1[i])){
        
        
        
        
        if (grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates_without_cluster1[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
          
        }
        
        
        if(grepl("^factor\\(", covariates_without_cluster1[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster1[i])
          covariate_value <- eval(parse(text =  covariate_name))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates_without_cluster1[i]) & !grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
          
          covariate_name <- covariates_without_cluster1[i]
          
          if(is.factor(eval(parse(text =  covariates_without_cluster1[i])))){
            covariate_value <- as.character(eval(parse(text =  covariates_without_cluster1[i])))
          }
          if(!is.factor(eval(parse(text =  covariates_without_cluster1[i])))){
            covariate_value <- eval(parse(text =  covariates_without_cluster1[i]))
          }
        }
        
        
        
        
        
        if(is.character(covariate_value)) {
          
          
          if(grepl("^factor\\(", covariates_without_cluster1[i])) {
            covfactors1<- c(covfactors1,gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster1[i]))
          }
          
          if(!grepl("^factor\\(", covariates_without_cluster1[i]) & !grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
            covfactors1<- c(covfactors1,covariates_without_cluster1[i])
          }
          
          column <- factor(covariate_value)
          column <- model.matrix(~ column -1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          column <- as.matrix(column[,-1])
          
          
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates_without_cluster1[i])
        }
        
        colnames(column) <- columnname 
        Xmat1 <- cbind(Xmat1,column)
        
      }
      
      
      
      
      
      if(!grepl("\\$",covariates_without_cluster1[i])){
        
        
        if (grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates_without_cluster1[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
          
          
        }
        
        
        if(grepl("^factor\\(", covariates_without_cluster1[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster1[i])
          covariate_value <- eval(parse(text = paste0("data$", covariate_name)))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates_without_cluster1[i]) & !grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
          
          covariate_name <- covariates_without_cluster1[i]
          
          if(is.factor(eval(parse(text = paste0("data$", covariates_without_cluster1[i]))))){
            covariate_value <- as.character(eval(parse(text = paste0("data$", covariates_without_cluster1[i]))))
          }
          if(!is.factor(eval(parse(text = paste0("data$", covariates_without_cluster1[i]))))){
            covariate_value <- eval(parse(text = paste0("data$", covariates_without_cluster1[i])))
          }
          
        }
        
        if(is.character(covariate_value)) {
          
          if(grepl("^factor\\(", covariates_without_cluster1[i])) {
            covfactors1<- c(covfactors1,gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster1[i]))
          }
          
          if(!grepl("^factor\\(", covariates_without_cluster1[i]) & !grepl("^I\\(.*\\)$", covariates_without_cluster1[i])) {
            covfactors1<- c(covfactors1,covariates_without_cluster1[i])
          }
          
          
          column <- factor(covariate_value)
          column <- model.matrix(~ column-1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          
          column <- as.matrix(column[,-1])
          
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates_without_cluster1[i])
          
        }
        
        colnames(column) <- columnname 
        
        Xmat1 <- cbind(Xmat1,column)
        
        
      }
      
    }
  }
  
  
  
  ###XMAT2
  covfactors2 <- c()
  terms_object2 <- terms(formula.terminalEvent)
  covariates2 <- attr(terms_object2, "term.labels")
  
  
  
  
  if(length(covariates2)>0){
    
    Xmat2 <- NULL
    
    
    for(i in 1:length(covariates2)){
      
      columnname <- c()
      
      
      if(grepl("\\$",covariates2[i])){
        
        
        
        if (grepl("^I\\(.*\\)$", covariates2[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates2[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
         
          
        }
        
        
        
        if(grepl("^factor\\(", covariates2[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates2[i])
          covariate_value <- eval(parse(text =  covariate_name))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates2[i])& !grepl("^I\\(.*\\)$", covariates2[i])) {
          
          covariate_name <- covariates2[i]
          
          if(is.factor(eval(parse(text =  covariates2[i])))){
            covariate_value <- as.character(eval(parse(text =  covariates2[i])))
          }
          if(!is.factor(eval(parse(text =  covariates2[i])))){
            covariate_value <- eval(parse(text =  covariates2[i]))
          }
        }
        
        
        
        
        
        if(is.character(covariate_value)) {
          
          
          if(grepl("^factor\\(", covariates2[i])) {
            covfactors2<- c(covfactors2,gsub("^factor\\((.*)\\)$", "\\1", covariates2[i]))
          }
          
          if(!grepl("^factor\\(", covariates2[i])& !grepl("^I\\(.*\\)$", covariates2[i])) {
            covfactors2<- c(covfactors2,covariates2[i])
          }
          
          column <- factor(covariate_value)
          column <- model.matrix(~ column -1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          column <- as.matrix(column[,-1])
          
          
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates2[i])
        }
        
        colnames(column) <- columnname 
        Xmat2 <- cbind(Xmat2,column)
        
      }
      
      
      
      
      
      if(!grepl("\\$",covariates2[i])){
        
        
        if (grepl("^I\\(.*\\)$", covariates2[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates2[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
          
          
        }
        
        
        if(grepl("^factor\\(", covariates2[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates2[i])
          covariate_value <- eval(parse(text = paste0("data$", covariate_name)))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates2[i]) & !grepl("^I\\(.*\\)$", covariates2[i])) {
          covariate_name <- covariates2[i]
          
          if(is.factor(eval(parse(text = paste0("data$", covariates2[i]))))){
            covariate_value <- as.character(eval(parse(text = paste0("data$", covariates2[i]))))
          }
          if(!is.factor(eval(parse(text = paste0("data$", covariates2[i]))))){
            covariate_value <- eval(parse(text = paste0("data$", covariates2[i])))
          }
          
        }
        
        if(is.character(covariate_value)) {
          
          if(grepl("^factor\\(", covariates2[i])) {
            covfactors2<- c(covfactors2,gsub("^factor\\((.*)\\)$", "\\1", covariates2[i]))
          }
          
          if(!grepl("^factor\\(", covariates2[i])& !grepl("^I\\(.*\\)$", covariates2[i])) {
            covfactors2<- c(covfactors2,covariates2[i])
          }
          
          
          column <- factor(covariate_value)
          
          column <- model.matrix(~ column-1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          
          column <- as.matrix(column[,-1])
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates2[i])
          
        }
        
        colnames(column) <- columnname 
        
        Xmat2 <- cbind(Xmat2,column)
        
        
      }
      
    }
  }
  
  
  
  
  ###Xmat3
  covfactors3 <- c()
  terms_object3 <- terms(formula.terminalEvent)
  covariates3 <- attr(terms_object3, "term.labels")
  
  
  
  
  if(length(covariates3)>0){
    
    Xmat3 <- NULL
    
    
    for(i in 1:length(covariates3)){
      
      columnname <- c()
      
      
      if(grepl("\\$",covariates3[i])){
        
        
        
        if (grepl("^I\\(.*\\)$", covariates3[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates3[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
         
          
        }
        
        
        
        if(grepl("^factor\\(", covariates3[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates3[i])
          covariate_value <- eval(parse(text =  covariate_name))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates3[i]) & !grepl("^I\\(.*\\)$", covariates3[i])) {
          
          covariate_name <- covariates3[i]
          
          if(is.factor(eval(parse(text =  covariates3[i])))){
            covariate_value <- as.character(eval(parse(text =  covariates3[i])))
          }
          if(!is.factor(eval(parse(text =  covariates3[i])))){
            covariate_value <- eval(parse(text =  covariates3[i]))
          }
        }
        
        
        
        
        
        if(is.character(covariate_value)) {
          
          
          if(grepl("^factor\\(", covariates3[i])) {
            covfactors3<- c(covfactors3,gsub("^factor\\((.*)\\)$", "\\1", covariates3[i]))
          }
          
          if(!grepl("^factor\\(", covariates3[i]) & !grepl("^I\\(.*\\)$", covariates3[i])) {
            covfactors3<- c(covfactors3,covariates3[i])
          }
          
          column <- factor(covariate_value)
          column <- model.matrix(~ column -1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          column <- as.matrix(column[,-1])
          
          
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates3[i])
        }
        
        colnames(column) <- columnname 
        Xmat3 <- cbind(Xmat3,column)
        
      }
      
      
      
      
      
      if(!grepl("\\$",covariates3[i])){
        
        
        if (grepl("^I\\(.*\\)$", covariates3[i])) {
          inner_expr <- gsub("^I\\((.*)\\)$", "\\1", covariates3[i])
          covariate_value <- eval(parse(text = inner_expr), envir = data)
          column <- as.matrix(covariate_value)
          
          
        }
        
        
        if(grepl("^factor\\(", covariates3[i])) {
          covariate_name <- gsub("^factor\\((.*)\\)$", "\\1", covariates3[i])
          covariate_value <- eval(parse(text = paste0("data$", covariate_name)))
          
          covariate_value <- as.character(covariate_value)
        }
        
        
        
        if(!grepl("^factor\\(", covariates3[i]) & !grepl("^I\\(.*\\)$", covariates3[i])) {
          covariate_name <- covariates3[i]
          
          if(is.factor(eval(parse(text = paste0("data$", covariates3[i]))))){
            covariate_value <- as.character(eval(parse(text = paste0("data$", covariates3[i]))))
          }
          if(!is.factor(eval(parse(text = paste0("data$", covariates3[i]))))){
            covariate_value <- eval(parse(text = paste0("data$", covariates3[i])))
          }
          
        }
        
        if(is.character(covariate_value)) {
          
          if(grepl("^factor\\(", covariates3[i])) {
            covfactors3<- c(covfactors3,gsub("^factor\\((.*)\\)$", "\\1", covariates3[i]))
          }
          
          if(!grepl("^factor\\(", covariates3[i]) & !grepl("^I\\(.*\\)$", covariates3[i])) {
            covfactors3<- c(covfactors3,covariates3[i])
          }
          
          
          column <- factor(covariate_value)
          column <- model.matrix(~ column-1)
          columnname <- paste0(covariate_name,paste0(".",gsub("column", "\\1", colnames(column))[-1]))
          column <- as.matrix(column[,-1])
          
          
        }
        
        if(!is.character(covariate_value)) {
          column <- as.matrix(covariate_value)
          columnname <- as.character(covariates3[i])
          
        }
        
        colnames(column) <- columnname 
        
        Xmat3 <- cbind(Xmat3,column)
        
        
      }
      
    }
  }
  

  
  
  
  
  
  terms_object <- terms(formula)
  covariates <- attr(terms_object, "term.labels")
  
  if(length(covariates)==1 && grepl("^cluster\\(.*\\)", covariates)){
    Xmat1 <- matrix(nrow=nrow(data),ncol=0)
  }
  
  
  if(formula[[3]]==1){
    Xmat1 <- matrix(nrow=nrow(data),ncol=0)
  } 
  
  
  if(formula.terminalEvent[[3]]==1){
    Xmat2 <- matrix(nrow=nrow(data),ncol=0)
    Xmat3 <- matrix(nrow=nrow(data),ncol=0)
  }
  
  
  dups <- duplicated(colnames(Xmat1))
  Xmat1 <- Xmat1[, !dups,drop=FALSE]
  covfactors1 <- unique(covfactors1)
  
  dups <- duplicated(colnames(Xmat2))
  Xmat2 <- Xmat2[, !dups,drop=FALSE]
  covfactors2 <- unique(covfactors2)
  
  dups <- duplicated(colnames(Xmat3))
  Xmat3 <- Xmat3[, !dups,drop=FALSE]
  covfactors3 <- unique(covfactors3)
  
  allcovariates <- c()
  
  if(length(colnames(Xmat1))>0){
    allcovariates <- c(allcovariates,paste0(colnames(Xmat1),paste0(".","01")))
  }
  
  if(length(colnames(Xmat2))>0){
    allcovariates <- c(allcovariates, paste0(colnames(Xmat2),paste0(".","02")))
  }
  
  
  if(length(colnames(Xmat3))>0){
    allcovariates <- c(allcovariates,paste0(colnames(Xmat3),paste0(".","12")))
  }
  
  
  
  
  
  ##
  if(length(covariates_without_cluster1)>0){
    covariate_terms1 <- paste0("Xmat1[,", seq_len(ncol(Xmat1)), "]", collapse = " + ")
    dynamic_formula1 <- as.formula(paste("Surv(l, y1, delta1) ~", covariate_terms1))
    
    
    suppressWarnings({frailtyinit1 <- frailtyPenal(dynamic_formula1,data=data_init,
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)
    })
  }
  
  
  if(length(covariates_without_cluster1)==0){
    
    suppressWarnings({frailtyinit1 <- frailtyPenal(Surv(l, y1, delta1)~1,data=data_init,
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)
    
    })
  }
  
  
  if(length(covariates2)>0){
    covariate_terms2 <- paste0("Xmat2[,", seq_len(ncol(Xmat2)), "]", collapse = " + ")
    dynamic_formula2 <- as.formula(paste("Surv(l, y2, delta2) ~", covariate_terms2))
    suppressWarnings({frailtyinit2 <- frailtyPenal(dynamic_formula2,data=data_init,
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)})
  }
  
  if(length(covariates2)==0){
    
    suppressWarnings({frailtyinit2 <- frailtyPenal(Surv(l, y2, delta2)~1,data=data_init,
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)})
  }
  
  
  sojourn <- y2[which(delta1==1)] - y1[which(delta1==1)]
  delta22 <- delta2[which(delta1==1)]
  
  
  
  if(length(covariates3)>0){
    Xmat33 <- as.matrix(Xmat3[which(delta1==1),])
    
    covariate_terms3 <- paste0("Xmat33[,", seq_len(ncol(Xmat33)), "]", collapse = " + ")
    dynamic_formula3 <- as.formula(paste("Surv(sojourn,delta22) ~", covariate_terms3))
    suppressWarnings({frailtyinit3 <- frailtyPenal(dynamic_formula3,data=data_init[which(delta1==1),],
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)})
  }
  
  
  
  if(length(covariates3)==0){
    
    suppressWarnings({frailtyinit3 <- frailtyPenal(Surv(sojourn,delta22)~1,data=data_init[which(delta1==1),],
                                                   hazard = "Weibull",maxit = maxit,print.times = FALSE)})
  }
  
  
  
  
  
  scale.weib1 <- exp(0.1)
  shape.weib1 <- exp(0.1)
  
  if(length(covariates_without_cluster1)>0){
    coef1 <- rep(0.1,ncol(Xmat1))
  }
  
  
  scale.weib2 <- exp(0.1)
  shape.weib2 <- exp(0.1)
  
  if(length(covariates2)>0){
    coef2 <- rep(0.1,ncol(Xmat2))
  }
  
  scale.weib3 <- exp(0.1)
  shape.weib3 <- exp(0.1)
  
  if(length(covariates3)>0){
    coef3 <- rep(0.1,ncol(Xmat3))
  }
  
  if(frailtyinit1$istop==1){
    scale.weib1 <- frailtyinit1$scale.weib[1]
    shape.weib1 <- frailtyinit1$shape.weib[1]
    
    if(length(covariates_without_cluster1)>0){
      coef1 <- frailtyinit1$coef
    }
  }
  
  
  
  
  
  if(frailtyinit2$istop==1){
    scale.weib2 <- frailtyinit2$scale.weib[1]
    shape.weib2 <- frailtyinit2$shape.weib[1]
    
    if(length(covariates2)>0){
      coef2 <- frailtyinit2$coef
    }
  }
  
  
  
  
  
  if(frailtyinit3$istop==1){
    scale.weib3 <- frailtyinit3$scale.weib[1]
    shape.weib3 <- frailtyinit3$shape.weib[1]
    
    if(length(covariates3)>0){
      coef3 <- frailtyinit3$coef
    }
    
  }
  
  
  
  
  startVals     <- c(log(scale.weib1), log(shape.weib1),
                     log(scale.weib2), log(shape.weib2),
                     log(scale.weib3), log(shape.weib3))
  
  hazard_coef <- c((scale.weib1), (shape.weib1),
                   (scale.weib2), (shape.weib2),
                   (scale.weib3), (shape.weib3))
  
  terms_object <- terms(formula)
  covariates <- attr(terms_object, "term.labels")
  
  
  if(any(grepl("^cluster\\(.*\\)", covariates))){ 
    
    
    startVals <- c(startVals, (sqrt(0.1)))
    
  }
  
  
  if(length(covariates_without_cluster1)>0){
    startVals <- c(startVals,as.numeric(coef1))
  }
  
  if(length(covariates2)>0){
    startVals <- c(startVals,as.numeric(coef2))
  }
  
  
  if(length(covariates3)>0){
    startVals <- c(startVals,as.numeric(coef3))
  }
  
  
  if(!missing(init.hazard.weib)){
    
    init.hazard.weib=init.hazard.weib
    
   
    startVals[1,3,5]<- log(c(as.numeric(unlist(init.hazard.weib[1,3,5]))))
    startVals[2,4,6]<- log(c(as.numeric(unlist(init.hazard.weib[2,4,6]))))
    
  }
  
  
  if(!missing(init.B)){
    
    init.B=init.B
    
    if(any(grepl("^cluster\\(.*\\)", covariates))){
      
      
      if(length(covariates_without_cluster1)>0){
        
        init.B1 <- init.B[1:ncol(Xmat1)]
        
       
        
        
        startVals[8:(8+ncol(Xmat1)-1)] <- c(as.numeric(unlist(init.B1)))
        
        
        
      }
      
      if(length(covariates2)>0){
        
        
        init.B2 <- init.B[(ncol(Xmat1)+1):(ncol(Xmat1)+ncol(Xmat2))]
        
       
        
        startVals[(8+ncol(Xmat1)):(8+ncol(Xmat1)+ncol(Xmat2)-1)]<- c(as.numeric(unlist(init.B2)))
        
      }
      
      
      if(length(covariates3)>0){
        
        
        init.B3 <- init.B[(ncol(Xmat1)+ncol(Xmat2)+1):(ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))]
        
            
        startVals[(8+ncol(Xmat1)+ncol(Xmat2)):(8+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]<- c(as.numeric(unlist(init.B3)))
        
      }
      
      
      
    }
    
    if(!any(grepl("^cluster\\(.*\\)", covariates))){
      
      
      if(length(covariates_without_cluster1)>0){
        
        init.B1 <- init.B[1:ncol(Xmat1)]
        
       
        
        startVals[7:(7+ncol(Xmat1)-1)]<- c(as.numeric(unlist(init.B1)))
        
        
        
      }
      
      if(length(covariates2)>0){
        
        
        init.B2 <- init.B[(ncol(Xmat1)+1):(ncol(Xmat1)+ncol(Xmat2))]
       
        
        startVals[(7+ncol(Xmat1)):(7+ncol(Xmat1)+ncol(Xmat2)-1)]<- c(as.numeric(unlist(init.B2)))
        
      }
      
      
      if(length(covariates3)>0){
        
        
        init.B3 <- init.B[(ncol(Xmat1)+ncol(Xmat2)+1):(ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))]
        
          
        startVals[(7+ncol(Xmat1)+ncol(Xmat2)):(7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]<- c(as.numeric(unlist(init.B3)))
        
      }
      
      
      
      
      
    }
    
  }
  
  if(missing(init.Theta)){
    init.Theta <- sqrt(0.1)
    startVals[5] <- init.Theta
  }
  
  if(any(grepl("^cluster\\(.*\\)", covariates))){ 
    
    if(!missing(init.Theta)){
      init.Theta=init.Theta
      startVals[7] <- init.Theta
    }
  }
  
  
  
  
  
  
  ##################
  
  if (!missing(partialH)) {
    partialH <- unique(partialH)
    
    if (!all(partialH == floor(partialH)) || any(partialH <= 0)) {
      stop("'partialH' should be a vector of positive integers.")
    }
    
    valid_indices <- 1:length(startVals)  
    
    
    if (any(!partialH %in% valid_indices)) {
      stop("'partialH' contains invalid indices. Only indices within 1 to ", length(startVals), " are allowed.")
    }
  }
  
  if(missing(partialH)){
    partialH <- c()
  }
  
  
  
  ##############
  
  if(model == "Semi-Markov")
  {
    
    
    
    
    if(any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
     
      cat("Be patient the program is computing...\n")
      logLike <- function(p) logLike.group.weibull.SCR.SM.LT.FRAILTY.SEMI.MARKOV(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2,l=l,
                                                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,data=data,group=group)
      error <- NULL
      tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                   ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                   epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
               error = function(e) {
                 error <<- e$message
               })
      
      
    }
    
    
    
    
    
    
    
    
    
    if(!any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
     
      cat("Be patient the program is computing...\n")
      logLike <- function(p) logLike.weibull.SCR.SM.LT.SANS.FRAILTY.SEMI.MARKOV(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3)
      
      error <- NULL
      
      tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                   ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                   epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
               error = function(e) {
                 error <<- e$message
               })
      
      
      
      
    }
    
  }
  ##
  
  
  if(model == "Markov")
  {
    
    
    
    
    if(any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
   
      cat("Be patient the program is computing...\n")
      logLike <- function(p) logLike.group.weibull.SCR.SM.LT.FRAILTY.MARKOV(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2,l=l,
                                                                            
                                                                            Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,data=data,group=group)
      
      error <- NULL
      tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                   ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                   epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
               error = function(e) {
                 error <<- e$message
               })
      
    }
    
    
    
    if(!any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
     
      cat("Be patient the program is computing...\n")
      logLike <- function(p) logLike.weibull.SCR.SM.LT.SANS.FRAILTY.MARKOV(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                           Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3)
      
      
      error <- NULL
      tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                   ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                   epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
               error = function(e) {
                 error <<- e$message
               })
      
    }
    
  }
  ##
  
  
  
  if (fit1$istop == 4) {
    stop("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
  } 
  
  
  if(is.null(fit1)){
    stop("An unexpected issue occurred during optimization. 'marqLevAlg' did not return a result.")
  }
  
  
  
  original_call <- match.call()
  
  
  
  
  
  
  
  if(direct==TRUE){
    cat("\nCall:\n")
    print(original_call)
    cat("\n")
  }
  
  
 
  
   
  
  
  v <- fit1$v
  n <- nrow(data)
  p <- length(fit1$b)
  vcov_matrix <- matrix(0, p, p)
  vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
  
  vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
  
  if(Frailty==1){
    vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                          exp(fit1$b[3]),exp(fit1$b[4]),
                          exp(fit1$b[5]),exp(fit1$b[6]),
                          2*fit1$b[7],rep(1,ncol(Xmat1)),
                          rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                        7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                        7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                          exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                          exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                          2*fit1$b[7],rep(1,ncol(Xmat1)),
                                                                                          rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                        7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                        7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
  }
  
  
  if(Frailty==0){
    vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                          exp(fit1$b[3]),exp(fit1$b[4]),
                          exp(fit1$b[5]),exp(fit1$b[6])
                          ,rep(1,ncol(Xmat1)),
                          rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                        6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                        6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                          exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                          exp(fit1$b[5]),exp(fit1$b[6])
                                                                                          ,rep(1,ncol(Xmat1)),
                                                                                          rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                        6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                        6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
  }
  
  
  
  
  
  
  scale.weib <- exp(fit1$b[c(1,3,5)])
  shape.weib <- exp(fit1$b[c(2,4,6)])
  vcov <- vcov_matrix
  call <- original_call
  crit <- fit1$istop
  grad <- fit1$grad
  if(Frailty==1){
    groups <- unique(group)
  }
  n.events <- c(length(which(delta1 == 1)),length(which(delta1 == 0 & delta2 == 1)),
                length(which(delta1 == 1 & delta2 == 1)),length(which(delta1 == 0 & delta2 == 0)))
  loglik <- fit1$fn.value
  
  if(Frailty==1){
    
    coef <- fit1$b[8:p]
    names(coef) <- allcovariates
  }
  if(Frailty==0){
    coef <- fit1$b[7:p]
    names(coef) <- allcovariates
  }
  n.iter <- fit1$ni
  AIC <- (1/n) *(length(fit1$b) - fit1$fn.value)
  
  
  
  beta_p.value <- c()
  terms_object <- terms(formula)
  covariates <- attr(terms_object, "term.labels")
  covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
  debut_reg_01 <- ifelse(Frailty == 1, 8, 7)
  if(length(covariates_without_cluster01)>0){
    
    covariates_transitions_01 <- fit1$b[debut_reg_01:(debut_reg_01 + ncol(Xmat1)-1)]
    
    for(i in 1:length(colnames(Xmat1))){
      beta_p.value<- c(beta_p.value,ifelse(covariates_transitions_01[i] > 0,
                                           2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                                           2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE))))
    }
    
    
  }
  
  terms_object <- terms(formula.terminalEvent)
  covariates <- attr(terms_object, "term.labels")
  covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
  debut_reg_02 <- debut_reg_01 + ncol(Xmat1)
  if(length(covariates_without_cluster02)>0){
    
    covariates_transitions_02 <- fit1$b[debut_reg_02:(debut_reg_02 + ncol(Xmat2)-1)]
    for(i in 1:length(colnames(Xmat2))){
      beta_p.value<- c(beta_p.value,ifelse(covariates_transitions_02[i] > 0,
                                           2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                                           2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE))))
    }
    
    
  }
  
  
  terms_object <- terms(formula.terminalEvent)
  covariates <- attr(terms_object, "term.labels")
  covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
  debut_reg_12 <- debut_reg_02 + ncol(Xmat2)
  if(length(covariates_without_cluster12)>0){
    
    covariates_transitions_12 <- fit1$b[debut_reg_12:(debut_reg_12 + ncol(Xmat3)-1)]
    for(i in 1:length(colnames(Xmat3))){
      beta_p.value<- c(beta_p.value,ifelse(covariates_transitions_12[i] > 0,
                                           2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                                           2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE))))
    }
    
  }
  
  if(Frailty==1){
    theta <- ((fit1$b[7])^2)
    VarTheta <- sqrt(vcov_matrix[7,7])
  }
  
  
  if(Frailty == 1){
    theta_p.value <- 1 - pnorm(theta/sqrt(VarTheta))
  }
  b <- c()
  b <- c(b,scale.weib[1],shape.weib[1],scale.weib[2],shape.weib[2],
         scale.weib[3],shape.weib[3])
  
  if(Frailty==1){
    b <- c(b,theta)
  }
  
  b <- c(b,coef)
  
  
  npar <- length(b)
  nvar=ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)
  
  
  ##BASELINE HAZARDS AND CONFIDENCE BANDS IF POSSIBLE
  
  if(1 %in% partialH | 2 %in% partialH){
    lam01 <- base_hazard(x01,scale.weib[1],shape.weib[1])
  }
  
  if(3 %in% partialH | 4 %in% partialH){
    lam02 <- base_hazard(x02,scale.weib[2],shape.weib[2])
  }
  
  if(5 %in% partialH | 6 %in% partialH){
    lam12 <- base_hazard(x12,scale.weib[3],shape.weib[3])
  }
  
  
  
  if(!(1 %in% partialH | 2 %in% partialH)){
    error_fit_haz_01 <- NULL
    tryCatch({
      
      lam01 <- hazard_and_confidence_bands_illdeath(x01, scale.weib[1], shape.weib[1], vcov_matrix[1:2, 1:2])
    }, error = function(e) {
      error_fit_haz_01 <<- e$message
    })
    
    if(!is.null(error_fit_haz_01)){
      lam01 <- base_hazard(x01, scale.weib[1], shape.weib[1])
    }
  }
  
  if(!(3 %in% partialH | 4 %in% partialH)){
    error_fit_haz_02 <- NULL
    tryCatch({
      lam02 <- hazard_and_confidence_bands_illdeath(x02, scale.weib[2], shape.weib[2], vcov_matrix[3:4, 3:4])
    }, error = function(e) {
      error_fit_haz_02 <<- e$message
    })
    
    if(!is.null(error_fit_haz_02)){
      lam02 <- base_hazard(x01, scale.weib[2], shape.weib[2])
    }
  }
  
  if(!(5 %in% partialH | 6 %in% partialH)){
    error_fit_haz_12 <- NULL
    tryCatch({
      lam12 <- hazard_and_confidence_bands_illdeath(x12, scale.weib[3], shape.weib[3], vcov_matrix[5:6, 5:6])
    }, error = function(e) {
      error_fit_haz_12 <<- e$message
    })
    
    if(!is.null(error_fit_haz_12)){
      lam12 <- base_hazard(x12, scale.weib[3], shape.weib[3])
    }
  }
  
  
  
  ####
  
  ##BASELINE SURVIVALS AND CONFIDENCE BANDS IF POSSIBLE
  
  if(1 %in% partialH | 2 %in% partialH){
    surv01 <- su(x01,scale.weib[1],shape.weib[1])
    
    median01 <- minmin(surv01[,2],x01)
    
    median.01 <- c(median.01=median01)
  }
  
  if(3 %in% partialH | 4 %in% partialH){
    surv02 <- su(x02,scale.weib[2],shape.weib[2])
    
    median02 <- minmin(surv02[,2],x02)
    
    median.02 <- c(median.02=median02)
  }
  
  if(5 %in% partialH | 6 %in% partialH){
    surv12 <- su(x12,scale.weib[3],shape.weib[3])
    
    median12 <- minmin(surv12[,2],x12)
    
    median.12 <- c(median.12=median12)
  }
  
  
  
  if(!(1 %in% partialH | 2 %in% partialH)){
    error_fit_surv_01 <- NULL
    tryCatch({
      
      surv01 <- survival_and_confidence_bands_illdeath(x01, scale.weib[1], shape.weib[1], vcov_matrix[1:2, 1:2])
    }, error = function(e) {
      error_fit_surv_01 <<- e$message
    })
    
    if(is.null(error_fit_surv_01)){
      
      median01 <- minmin(surv01[,2],x01)
      
      lower01 <- minmin(surv01[,3],x01)
      
      upper01 <- minmin(surv01[,4],x01)
      
      median.01 <- c(lower.01=lower01,median.01=median01,upper.01=upper01)
    }
    
    
    if(!is.null(error_fit_surv_01)){
      surv01 <- su(x01, scale.weib[1], shape.weib[1])
      
      median01 <- minmin(surv01[,2],x01)
      
      median.01 <- c(median.01=median01)
      
      
    }
  }
  
  if(!(3 %in% partialH | 4 %in% partialH)){
    error_fit_surv_02 <- NULL
    tryCatch({
      
      surv02 <- survival_and_confidence_bands_illdeath(x02, scale.weib[2], shape.weib[2], vcov_matrix[3:4, 3:4])
    }, error = function(e) {
      error_fit_surv_02 <<- e$message
    })
    
    
    if(is.null(error_fit_surv_02)){
      
      median02 <- minmin(surv02[,2],x02)
      
      lower02 <- minmin(surv02[,3],x02)
      
      upper02 <- minmin(surv02[,4],x02)
      
      median.02 <- c(lower.02=lower02,median.02=median02,upper.02=upper02)
    }
    
    if(!is.null(error_fit_surv_02)){
      surv02 <- su(x02, scale.weib[2], shape.weib[2])
      
      median02 <- minmin(surv02[,2],x02)
      
      median.02 <- c(median.02=median02)
    } 
    
  }
  
  if(!(5 %in% partialH | 6 %in% partialH)){
    error_fit_surv_12 <- NULL
    tryCatch({
      
      surv12 <- survival_and_confidence_bands_illdeath(x12, scale.weib[3], shape.weib[3], vcov_matrix[5:6, 5:6])
    }, error = function(e) {
      error_fit_surv_12 <<- e$message
    })
    
    
    if(is.null(error_fit_surv_12)){
      
      median12 <- minmin(surv12[,2],x12)
      
      lower12 <- minmin(surv12[,3],x12)
      
      upper12 <- minmin(surv12[,4],x12)
      
      median.12 <- c(lower.12=lower12,median.12=median12,upper.12=upper12)
    }
    
    if(!is.null(error_fit_surv_12)){
      surv12 <- su(x12, scale.weib[3], shape.weib[3])
      
      median12 <- minmin(surv12[,2],x12)
      
      median.12 <- c(median.12=median12)
    }  
    
  }
  
  ###
  
  ############### FRAILTY PREDICTION
  
  if(Frailty==1){
    
    eta1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% covariates_transitions_01)
    eta2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% covariates_transitions_02)
    eta3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% covariates_transitions_12)
    
    
    q.y1 <- exp(eta1)*(y1/scale.weib[1])^shape.weib[1]
    q.y2 <- exp(eta2)*(y2/scale.weib[2])^shape.weib[2]
    
    if (model == "Semi-Markov") {
      q.y3 <- ifelse(delta1 == 1, exp(eta3) * ((y2 - y1) / scale.weib[3])^shape.weib[3], 0)
    }
    
    
    if (model == "Markov") {
      q.y3 <- ifelse(delta1 == 1, exp(eta3) * ((y2) / scale.weib[3])^shape.weib[3], 0)
    }
    
    if(length(surv_trunc1)==3){
      q.y1 <- q.y1 - exp(eta1)*(l/scale.weib[1])^shape.weib[1]
      q.y2 <- q.y2 - exp(eta2)*(l/scale.weib[2])^shape.weib[2]
    }
    
    
    
    
    frailty.pred <- c() 
    frailty.var <- c() 
    frailty.sd <- c() 
    
    
    
    for(j in unique(data$group)){
      cluster <- which(data$group==j)
      
      
      num.pred <- sum(delta1[cluster])+sum(delta2[cluster])+(1/theta)
      denom.pred <- sum(q.y1[cluster])+ sum(q.y2[cluster])+ sum(q.y3[cluster])+(1/theta)
      frailty.pred <- c(frailty.pred,num.pred/denom.pred)
      frailty.var <- c(frailty.var,num.pred/(denom.pred)^2)
      frailty.sd <- c(frailty.sd,sqrt(frailty.var))
      
    }
    
  }
  
  
  linear.pred01 <- c()
  linear.pred02 <- c()
  linear.pred12 <- c()
  
  if(Frailty==1){
    
    for(i in 1:nrow(data)){
      
      cluster <- data$group[i]
      
      linear.pred01 <- c(linear.pred01,eta1[i]+log(frailty.pred[cluster]))
      linear.pred02 <- c(linear.pred02,eta2[i]+log(frailty.pred[cluster]))
      linear.pred12 <- c(linear.pred12,eta3[i]+log(frailty.pred[cluster]))
    }
  }
  
  
  if(Frailty==0){
    
    
    eta1 <- if (ncol(Xmat1) == 0) rep(0, nrow(Xmat1)) else as.vector(Xmat1 %*% covariates_transitions_01)
    eta2 <- if (ncol(Xmat2) == 0) rep(0, nrow(Xmat2)) else as.vector(Xmat2 %*% covariates_transitions_02)
    eta3 <- if (ncol(Xmat3) == 0) rep(0, nrow(Xmat3)) else as.vector(Xmat3 %*% covariates_transitions_12)
    
    
    
    
    
    
    for(i in 1:nrow(data)){
      
      
      
      linear.pred01 <- c(linear.pred01,eta1[i])
      linear.pred02 <- c(linear.pred02,eta2[i])
      linear.pred12 <- c(linear.pred12,eta3[i])
    }
  }
  
  
  

  
  
  ca= fit1$ca
  cb=fit1$cb
  rdm=fit1$rdm
  
  
  
  covfactors <- c(covfactors1,covfactors2,covfactors3)
  
  
  
  
  
  covfactors1 <- unique(covfactors1)
  covfactors2 <- unique(covfactors2)
  covfactors3 <- unique(covfactors3)
  covariates_without_cluster_without_factor_1 <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster1))
  covariates_without_factor_2 <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates2))
  covariates_without_factor_3 <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates3))
  
  
  
  
  if(length(covfactors1)>0){
    
    if(Frailty==TRUE){
      vcov_coef01 <- vcov_matrix[8:(8+ncol(Xmat1)-1),8:(8+ncol(Xmat1)-1)]
      coef01 <- b[8:(8+ncol(Xmat1)-1)]
    }
    
    if(Frailty==FALSE){
      vcov_coef01 <- vcov_matrix[7:(7+ncol(Xmat1)-1),7:(7+ncol(Xmat1)-1)]
      coef01 <- b[7:(7+ncol(Xmat1)-1)]
    }
    
    global_chisq.01 <- c()
    dof_chisq.01 <- c()
    global_chisq.test.01 <- c()
    p.global_chisq.01 <- c()
    
    for(factor in covfactors1){
      idx_factor <- which(covariates_without_cluster_without_factor_1==factor)
      
      if(grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text =  factor))))) -1
      }
      
      if(!grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text = paste0("data$", factor))))))-1
      }
      
      
      
      if(ncol(Xmat1)>1){
        W <- t(coef01[(idx_factor):(idx_factor+indic_factor-1)])%*% solve(vcov_coef01[(idx_factor):(idx_factor+indic_factor-1),(idx_factor):(idx_factor+indic_factor-1)])%*%
          coef01[(idx_factor):(idx_factor+indic_factor-1)]
      }
      
      if(ncol(Xmat1)==1){
        W <- (coef01^2)/(vcov_coef01)
      }
      
      p <- 1-pchisq(W,df=indic_factor)
      
      
      
      
      
      global_chisq.01 <- c(global_chisq.01,W)
      dof_chisq.01 <- c(dof_chisq.01,indic_factor)
      p.global_chisq.01 <- c(p.global_chisq.01,p)
    }
    
    global_chisq.test.01 <- ifelse(length(global_chisq.01)>0,1,0)
    
    names(global_chisq.01) <- covfactors1
    
    names(dof_chisq.01) <- covfactors1
    
    names(p.global_chisq.01) <- covfactors1
    
  }
  
  if(length(covfactors1)==0){
    global_chisq.01 <- NULL
    dof_chisq.01 <- NULL
    p.global_chisq.01 <- NULL
    global_chisq.test.01 <- 0
  }
  
  
  
  
  if(length(covfactors2)>0){
    
    if(Frailty==TRUE){
      vcov_coef02 <- vcov_matrix[(8+ncol(Xmat1)):(8+ncol(Xmat1)+ncol(Xmat2)-1),(8+ncol(Xmat1)):(8+ncol(Xmat1)+ncol(Xmat2)-1)]
      coef02 <- b[(8+ncol(Xmat1)):(8+ncol(Xmat1)+ncol(Xmat2)-1)]
    }
    
    if(Frailty==FALSE){
      vcov_coef02 <- vcov_matrix[(7+ncol(Xmat1)):(7+ncol(Xmat1)+ncol(Xmat2)-1),(7+ncol(Xmat1)):(7+ncol(Xmat1)+ncol(Xmat2)-1)]
      coef02 <- b[(7+ncol(Xmat1)):(7+ncol(Xmat1)+ncol(Xmat2)-1)]
    }
    
    global_chisq.02 <- c()
    dof_chisq.02 <- c()
    global_chisq.test.02 <- c()
    p.global_chisq.02 <- c()
    
    for(factor in covfactors2){
      idx_factor <- which(covariates_without_factor_2==factor)
      
      if(grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text =  factor))))) -1
      }
      
      if(!grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text = paste0("data$", factor))))))-1
      }
      
      
      if(ncol(Xmat2)>1){
        W <- t(coef02[(idx_factor):(idx_factor+indic_factor-1)])%*% solve(vcov_coef02[(idx_factor):(idx_factor+indic_factor-1),(idx_factor):(idx_factor+indic_factor-1)])%*%
          coef02[(idx_factor):(idx_factor+indic_factor-1)]
      }
      
      if(ncol(Xmat2)==1){
        W <- (coef02^2)/(vcov_coef02)
      }
      p <- 1-pchisq(W,df=indic_factor)
      
      
      
      
      
      global_chisq.02 <- c(global_chisq.02,W)
      dof_chisq.02 <- c(dof_chisq.02,indic_factor)
      p.global_chisq.02 <- c(p.global_chisq.02,p)
    }
    
    global_chisq.test.02 <- ifelse(length(global_chisq.02)>0,1,0)
    
    names(global_chisq.02) <- covfactors2
    
    names(dof_chisq.02) <- covfactors2
    
    names(p.global_chisq.02) <- covfactors2
    
  }
  
  
  if(length(covfactors2)==0){
    global_chisq.02 <- NULL
    dof_chisq.02 <- NULL
    p.global_chisq.02 <- NULL
    global_chisq.test.02 <- 0
  }
  
  
  if(length(covfactors3)>0){
    
    if(Frailty==TRUE){
      vcov_coef12 <- vcov_matrix[(8+ncol(Xmat1)+ncol(Xmat2)):(8+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1),(8+ncol(Xmat1)+ncol(Xmat2)):(8+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]
      coef12 <- b[(8+ncol(Xmat1)+ncol(Xmat2)):(8+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]
    }
    
    if(Frailty==FALSE){
      vcov_coef12 <- vcov_matrix[(7+ncol(Xmat1)+ncol(Xmat2)):(7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1),(7+ncol(Xmat1)+ncol(Xmat2)):(7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]
      coef12 <- b[(7+ncol(Xmat1)+ncol(Xmat2)):(7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)-1)]
    }
    
    global_chisq.12 <- c()
    dof_chisq.12 <- c()
    global_chisq.test.12 <- c()
    p.global_chisq.12 <- c()
    
    for(factor in covfactors3){
      idx_factor <- which(covariates_without_factor_3==factor)
      
      if(grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text =  factor))))) -1
      }
      
      if(!grepl("\\$",factor)){
        indic_factor <- length(levels(as.factor(eval(parse(text = paste0("data$", factor))))))-1
      }
      
      
      if(ncol(Xmat3)>1){
        W <- t(coef12[(idx_factor):(idx_factor+indic_factor-1)])%*% solve(vcov_coef12[(idx_factor):(idx_factor+indic_factor-1),(idx_factor):(idx_factor+indic_factor-1)])%*%
          coef12[(idx_factor):(idx_factor+indic_factor-1)]
      }
      
      if(ncol(Xmat3)==1){
        W <- (coef12^2)/(vcov_coef12)
      }
      
      p <- 1-pchisq(W,df=indic_factor)
      
      
      
      
      
      global_chisq.12 <- c(global_chisq.12,W)
      dof_chisq.12 <- c(dof_chisq.12,indic_factor)
      p.global_chisq.12 <- c(p.global_chisq.12,p)
    }
    
    global_chisq.test.12 <- ifelse(length(global_chisq.12)>0,1,0)
    
    names(global_chisq.12) <- covfactors3
    
    names(dof_chisq.12) <- covfactors3
    
    names(p.global_chisq.12) <- covfactors3
    
  }
  
  
  
  
  if(length(covfactors3)==0){
    global_chisq.12 <- NULL
    dof_chisq.12 <- NULL
    p.global_chisq.12 <- NULL
    global_chisq.test.12 <- 0
  }
  
  
  
  
  
  if(model == "Semi-Markov")
  {
    
    
    
    
    if(any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
      cat("Illness-death model with shared frailty between transitions\n")
      cat("Using Weibull baseline hazard functions \n")
      
      
      if( length(surv_terms1)==3){
        cat("Left truncation structure is used\n")
      }
      
      cat("Semi-Markov model is used for the transition 1->2\n")
      cat("\n")
    }
    
    if(!any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
      cat("Illness-death model\n")
      cat("Using Weibull baseline hazard functions \n")
      
      
      if( length(surv_terms1)==3){
        cat("Left truncation structure is used\n")
      }
      
      cat("Semi-Markov model is used for the transition 1->2\n")
      cat("\n")
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(model == "Markov")
  {
    
    
    
    
    if(any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
      cat("Illness-death model with shared frailty between transitions\n")
      cat("Using Weibull baseline hazard functions \n")
      
      
      if(length(surv_terms1)==3){
        cat("Left truncation structure is used\n")
      }
      
      cat("Markov model is used for the transition 1->2\n")
      cat("\n")
    }
    
    
    if(!any(grepl("cluster\\(.*\\)", deparse(formula)))){
      
      cat("Illness-death model\n")
      cat("Using Weibull baseline hazard functions \n")
      
      
      if(length(surv_terms1)==3){
        cat("Left truncation structure is used\n")
      }
      
      cat("Markov model is used for the transition 1->2\n")
      cat("\n")
    }
  }
  
  
  
  
  trunc <- ifelse(length(surv_terms1)==3,TRUE,FALSE)
  
  
  
  if(direct==TRUE){
    if(fit1$istop==1){
      v <- fit1$v
      n <- length(fit1$b)
      vcov_matrix <- matrix(0, n, n)
      vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
      
      vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
      
      ### DELTA METHOD
      if(Frailty==1){
      vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                            exp(fit1$b[3]),exp(fit1$b[4]),
                            exp(fit1$b[5]),exp(fit1$b[6]),
                          2*fit1$b[7],rep(1,ncol(Xmat1)),
                          rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                          7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                          7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                            exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                            exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                            2*fit1$b[7],rep(1,ncol(Xmat1)),
                                                                                            rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                          7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                          7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
      }
      
      
      if(Frailty==0){
        vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                              exp(fit1$b[3]),exp(fit1$b[4]),
                              exp(fit1$b[5]),exp(fit1$b[6]),
                              rep(1,ncol(Xmat1)),
                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                              exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                              exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                              rep(1,ncol(Xmat1)),
                                                                                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
      }
      terms_object <- terms(formula)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_01 <- ifelse(Frailty == 1, 8, 7)
      
      if(length(covariates_without_cluster01)>0){
        covariates_transitions_01 <- fit1$b[debut_reg_01:(debut_reg_01 + ncol(Xmat1)-1)]
        
        cat("Transition 0 -> 1:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat1)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat1))){
          
          cov_name <- colnames(Xmat1)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_01[i], 
                      exp(covariates_transitions_01[i]),
                      sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), 
                      covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]),
                      ifelse(covariates_transitions_01[i] > 0,
                             2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        if (length(global_chisq.01) > 0) {
          # Filter factors to only include those with dof > 1, as per standard global test definition
          factors_with_mult_dof <- names(dof_chisq.01)[dof_chisq.01 > 1]
          
          if (length(factors_with_mult_dof) > 0) { # Only print header if there are such factors
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.01[factor_name],
                          dof_chisq.01[factor_name],
                          p.global_chisq.01[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      
      
      
      terms_object <- terms(formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_02 <- debut_reg_01 + ncol(Xmat1)
      if(length(covariates_without_cluster02)>0){
        
        covariates_transitions_02 <- fit1$b[debut_reg_02:(debut_reg_02 + ncol(Xmat2)-1)]
        
        cat("Transition 0 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat2)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat2))){
          
          cov_name <- colnames(Xmat2)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_02[i], 
                      exp(covariates_transitions_02[i]),
                      sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), 
                      covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]),
                      ifelse(covariates_transitions_02[i] > 0,
                             2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        if (length(global_chisq.02) > 0) {
          
          factors_with_mult_dof <- names(dof_chisq.02)[dof_chisq.02 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.02[factor_name],
                          dof_chisq.02[factor_name],
                          p.global_chisq.02[factor_name]
              ))
            }
            cat("\n") 
          }
        }
      }
      
      
      
      
      
      
      terms_object <- terms(formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_12 <- debut_reg_02 + ncol(Xmat2)
      if(length(covariates_without_cluster12)>0){
        
        covariates_transitions_12 <- fit1$b[debut_reg_12:(debut_reg_12 + ncol(Xmat3)-1)]
        
        cat("Transition 1 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat3)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat3))){
          
          cov_name <- colnames(Xmat3)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_12[i], 
                      exp(covariates_transitions_12[i]),
                      sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), 
                      covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]),
                      ifelse(covariates_transitions_12[i] > 0,
                             2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        if (length(global_chisq.12) > 0) {
          
          factors_with_mult_dof <- names(dof_chisq.12)[dof_chisq.12 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.12[factor_name],
                          dof_chisq.12[factor_name],
                          p.global_chisq.12[factor_name]
              ))
            }
            cat("\n") 
          }
        }
      }
      
      
      
      if(Frailty == 1){
        cat("Frailty Parameters:\n")
        cat("------------\n")
        cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", ((fit1$b[7])^2), sqrt(vcov_matrix[7,7]), 
                    1- pnorm((((fit1$b[7])^2)) / sqrt(vcov_matrix[7,7]))))
        cat("\n")
      }
      
      cat("Scales and shapes of the Weibull baseline hazard\n")
      cat("------------\n")
      cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
      
     
      cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[1]), sqrt(vcov_matrix[1,1]),
                  exp(fit1$b[2]), sqrt(vcov_matrix[2,2])))
      cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[3]), sqrt(vcov_matrix[3,3]),
                  exp(fit1$b[4]), sqrt(vcov_matrix[4,4])))
      cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[5]), sqrt(vcov_matrix[5,5]),
                  exp(fit1$b[6]), sqrt(vcov_matrix[6,6])))
      
      cat("\n")
     
      cat("The expression of the Weibull hazard function is: \n")
      cat("           'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' \n" )
      cat("The expression of the Weibull survival function is: \n")
      cat("           'S(t) = exp[- (t/scale)^shape]' \n" )
      
      cat("\n")
      
      cat("Marginal log-likelihood = ", fit1$fn.value, "\n")
      cat("AIC = Akaike information Criterion = ", (1/nrow(data)) * (length(fit1$b) - fit1$fn.value), "\n")
      cat("           'AIC = (1/n)[np - l(.)]' \n")
      
      cat("\n")
      
      cat("Number of subjects = ", nrow(data), "\n")
      cat("             0 -> 1 = ", length(which(delta1 == 1)), "\n")
      cat("             0 -> 2 = ", length(which(delta1 == 0 & delta2 == 1)), "\n")
      cat("             1 -> 2 = ", length(which(delta1 == 1 & delta2 == 1)), "\n")
      cat("Lost to follow-up  = ", length(which(delta1 == 0 & delta2 == 0)), "\n")
      cat("\n")
      
      cat("Number of iterations: ", fit1$ni, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
      cat("All criteria were satisfied")
      
    }
  }
  
  if(direct==TRUE){
    if(fit1$istop==2){
      v <- fit1$v
      n <- length(fit1$b)
      vcov_matrix <- matrix(0, n, n)
      vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
      
      vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
      
      
      
      
      
      if(Frailty==1){
        vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                              exp(fit1$b[3]),exp(fit1$b[4]),
                              exp(fit1$b[5]),exp(fit1$b[6]),
                              2*fit1$b[7],rep(1,ncol(Xmat1)),
                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                            7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                            7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                              exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                              exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                              2*fit1$b[7],rep(1,ncol(Xmat1)),
                                                                                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                            7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                            7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
      }
      
      
      if(Frailty==0){
        vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                              exp(fit1$b[3]),exp(fit1$b[4]),
                              exp(fit1$b[5]),exp(fit1$b[6]),
                              rep(1,ncol(Xmat1)),
                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                              exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                              exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                              rep(1,ncol(Xmat1)),
                                                                                              rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                            6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
      }
      
      
      
      terms_object <- terms(formula)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_01 <- ifelse(Frailty == 1, 8, 7)
      if(length(covariates_without_cluster01)>0){
        
        covariates_transitions_01 <- fit1$b[debut_reg_01:(debut_reg_01 + ncol(Xmat1)-1)]
        
        cat("Transition 0 -> 1:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat1)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat1))){
          
          cov_name <- colnames(Xmat1)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_01[i], 
                      exp(covariates_transitions_01[i]),
                      sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), 
                      covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]),
                      ifelse(covariates_transitions_01[i] > 0,
                             2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        if (length(global_chisq.01) > 0) {
          
          factors_with_mult_dof <- names(dof_chisq.01)[dof_chisq.01 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.01[factor_name],
                          dof_chisq.01[factor_name],
                          p.global_chisq.01[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      
      
      terms_object <- terms(formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_02 <- debut_reg_01 + ncol(Xmat1)
      if(length(covariates_without_cluster02)>0){
        
        covariates_transitions_02 <- fit1$b[debut_reg_02:(debut_reg_02 + ncol(Xmat2)-1)]
        
        cat("Transition 0 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat2)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat2))){
          
          cov_name <- colnames(Xmat2)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_02[i], 
                      exp(covariates_transitions_02[i]),
                      sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), 
                      covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]),
                      ifelse(covariates_transitions_02[i] > 0,
                             2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        
        if (length(global_chisq.02) > 0) {
          
          factors_with_mult_dof <- names(dof_chisq.02)[dof_chisq.02 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.02[factor_name],
                          dof_chisq.02[factor_name],
                          p.global_chisq.02[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      
      
      terms_object <- terms(formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_12 <- debut_reg_02 + ncol(Xmat2)
      if(length(covariates_without_cluster12)>0){
        covariates_transitions_12 <- fit1$b[debut_reg_12:(debut_reg_12 + ncol(Xmat3)-1)]
        
        cat("Transition 1 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(Xmat3)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
       
        for(i in 1:length(colnames(Xmat3))){
          
          cov_name <- colnames(Xmat3)[i]
          
          
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                      max_cov_name_length, cov_name, 
                      covariates_transitions_12[i], 
                      exp(covariates_transitions_12[i]),
                      sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), 
                      covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]),
                      ifelse(covariates_transitions_12[i] > 0,
                             2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                             2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
        }
        
        cat("\n")
        
        if (length(global_chisq.12) > 0) {
          
          factors_with_mult_dof <- names(dof_chisq.12)[dof_chisq.12 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq.12[factor_name],
                          dof_chisq.12[factor_name],
                          p.global_chisq.12[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      # Frailty Parameters
      if(Frailty == 1){
        cat("Frailty Parameters:\n")
        cat("------------\n")
        cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", ((fit1$b[7])^2), sqrt(vcov_matrix[7,7]), 
                    1- pnorm((((fit1$b[7])^2)) / sqrt(vcov_matrix[7,7]))))
        cat("\n")
      }
      
      # Weibull baseline hazard parameters
      cat("Scales and shapes of the Weibull baseline hazard\n")
      cat("------------\n")
      cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
      
     
      cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[1]), sqrt(vcov_matrix[1,1]),
                  exp(fit1$b[2]), sqrt(vcov_matrix[2,2])))
      cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[3]), sqrt(vcov_matrix[3,3]),
                  exp(fit1$b[4]), sqrt(vcov_matrix[4,4])))
      cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[5]), sqrt(vcov_matrix[5,5]),
                  exp(fit1$b[6]), sqrt(vcov_matrix[6,6])))
      
      cat("\n")
      
      cat("The expression of the Weibull hazard function is: \n")
      cat("           'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' \n" )
      cat("The expression of the Weibull survival function is: \n")
      cat("           'S(t) = exp[- (t/scale)^shape]' \n" )
      
      cat("\n")
      
      # Marginal log-likelihood and AIC
      cat("Marginal log-likelihood = ", fit1$fn.value, "\n")
      cat("AIC = Akaike information Criterion = ", (1/nrow(data)) * (length(fit1$b) - fit1$fn.value), "\n")
      cat("           'AIC = (1/n)[np - l(.)]' \n")
      
      cat("\n")
      
      # Number of subjects and transitions
      cat("Number of subjects = ", nrow(data), "\n")
      cat("             0 -> 1 = ", length(which(delta1 == 1)), "\n")
      cat("             0 -> 2 = ", length(which(delta1 == 0 & delta2 == 1)), "\n")
      cat("             1 -> 2 = ", length(which(delta1 == 1 & delta2 == 1)), "\n")
      cat("Lost to follow-up  = ", length(which(delta1 == 0 & delta2 == 0)), "\n")
      cat("\n")
      
      cat("Number of iterations: ", fit1$ni, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
      cat("Maximum number of iterations reached\n")
      if (fit1$ca > LIMparam) {
        cat("  Parameter criterion not reached. Threshold: ", LIMparam, "\n")
      }
      
      if (fit1$cb > LIMlogl) {
        cat("  Function criterion not reached. Threshold: ", LIMlogl, "\n")
      }
      
      if (fit1$rdm > LIMderiv) {
        cat("  Relative distance to optimum criterion not reached. Threshold: ", LIMderiv, "\n")
      }
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  if(direct==TRUE){
    if(length(partialH)>0){
      if(fit1$istop==3){
        v <- fit1$v
        n <- length(fit1$b)
        vcov_matrix <- matrix(0, n, n)
        vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
        
        vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
        
        
        if(Frailty==1){
          vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                exp(fit1$b[3]),exp(fit1$b[4]),
                                exp(fit1$b[5]),exp(fit1$b[6]),
                                2*fit1$b[7],rep(1,ncol(Xmat1)),
                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                                exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                                exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                                2*fit1$b[7],rep(1,ncol(Xmat1)),
                                                                                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
        }
        
        
        if(Frailty==0){
          vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                exp(fit1$b[3]),exp(fit1$b[4]),
                                exp(fit1$b[5]),exp(fit1$b[6]),
                                rep(1,ncol(Xmat1)),
                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                                exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                                exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                                rep(1,ncol(Xmat1)),
                                                                                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
        }
        
        
        terms_object <- terms(formula)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_01 <- ifelse(Frailty == 1, 8, 7)
        if(length(covariates_without_cluster01)>0){
          
          covariates_transitions_01 <- fit1$b[debut_reg_01:(debut_reg_01 + ncol(Xmat1)-1)]
          
          cat("Transition 0 -> 1:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat1)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat1))){
            
            cov_name <- colnames(Xmat1)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_01[i], 
                        exp(covariates_transitions_01[i]),
                        sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), 
                        covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]),
                        ifelse(covariates_transitions_01[i] > 0,
                               2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          if (length(global_chisq.01) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.01)[dof_chisq.01 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.01[factor_name],
                            dof_chisq.01[factor_name],
                            p.global_chisq.01[factor_name]
                ))
              }
              cat("\n") 
            }
          }
          
        }
        
        
        
        
        
        
        terms_object <- terms(formula.terminalEvent)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_02 <- debut_reg_01 + ncol(Xmat1)
        if(length(covariates_without_cluster02)>0){
          
          covariates_transitions_02 <- fit1$b[debut_reg_02:(debut_reg_02 + ncol(Xmat2)-1)]
          
          cat("Transition 0 -> 2:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat2)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat2))){
            
            cov_name <- colnames(Xmat2)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_02[i], 
                        exp(covariates_transitions_02[i]),
                        sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), 
                        covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]),
                        ifelse(covariates_transitions_02[i] > 0,
                               2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          
          if (length(global_chisq.02) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.02)[dof_chisq.02 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.02[factor_name],
                            dof_chisq.02[factor_name],
                            p.global_chisq.02[factor_name]
                ))
              }
              cat("\n") 
            }
          }
          
          
        }
        
        
        
        
        
        terms_object <- terms(formula.terminalEvent)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_12 <- debut_reg_02 + ncol(Xmat2)
        if(length(covariates_without_cluster12)>0){
          
          covariates_transitions_12 <- fit1$b[debut_reg_12:(debut_reg_12 + ncol(Xmat3)-1)]
          
          cat("Transition 1 -> 2:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat3)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat3))){
            
            cov_name <- colnames(Xmat3)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_12[i], 
                        exp(covariates_transitions_12[i]),
                        sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), 
                        covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]),
                        ifelse(covariates_transitions_12[i] > 0,
                               2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          
          if (length(global_chisq.12) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.12)[dof_chisq.12 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.12[factor_name],
                            dof_chisq.12[factor_name],
                            p.global_chisq.12[factor_name]
                ))
              }
              cat("\n") 
            }
          }
          
        }
        
        
        # Frailty Parameters
        if(Frailty == 1){
          cat("Frailty Parameters:\n")
          cat("------------\n")
          cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", ((fit1$b[7])^2), sqrt(vcov_matrix[7,7]), 
                      1- pnorm((((fit1$b[7])^2)) / sqrt(vcov_matrix[7,7]))))
          cat("\n")
        }
        
        # Weibull baseline hazard parameters
        cat("Scales and shapes of the Weibull baseline hazard\n")
        cat("------------\n")
        cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
        
       
        cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[1]), sqrt(vcov_matrix[1,1]),
                    exp(fit1$b[2]), sqrt(vcov_matrix[2,2])))
        cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[3]), sqrt(vcov_matrix[3,3]),
                    exp(fit1$b[4]), sqrt(vcov_matrix[4,4])))
        cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[5]), sqrt(vcov_matrix[5,5]),
                    exp(fit1$b[6]), sqrt(vcov_matrix[6,6])))
        
        cat("\n")
        
        cat("The expression of the Weibull hazard function is: \n")
        cat("           'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' \n" )
        cat("The expression of the Weibull survival function is: \n")
        cat("           'S(t) = exp[- (t/scale)^shape]' \n" )
        
        cat("\n")
        
        # Marginal log-likelihood and AIC
        cat("Marginal log-likelihood = ", fit1$fn.value, "\n")
        cat("AIC = Akaike information Criterion = ", (1/nrow(data)) * (length(fit1$b) - fit1$fn.value), "\n")
        cat("           'AIC = (1/n)[np - l(.)]' \n")
        
        cat("\n")
        
        # Number of subjects and transitions
        cat("Number of subjects = ", nrow(data), "\n")
        cat("             0 -> 1 = ", length(which(delta1 == 1)), "\n")
        cat("             0 -> 2 = ", length(which(delta1 == 0 & delta2 == 1)), "\n")
        cat("             1 -> 2 = ", length(which(delta1 == 1 & delta2 == 1)), "\n")
        cat("Lost to follow-up  = ", length(which(delta1 == 0 & delta2 == 0)), "\n")
        cat("\n")
        
        cat("Number of iterations: ", fit1$ni, "\n")
        cat("Convergence criteria:\n")
        cat("------------\n")
        cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
        cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
        cat("All criteria were satisfied, parameters in partialH were dropped from Hessian to define the relative distance to optimum")
        
      }
      
    }
  }
  
  if(direct==TRUE){
    if(length(partialH)>0){
      if(fit1$istop==2){
        v <- fit1$v
        n <- length(fit1$b)
        vcov_matrix <- matrix(0, n, n)
        vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
        
        vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
        
        
        
        if(Frailty==1){
          vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                exp(fit1$b[3]),exp(fit1$b[4]),
                                exp(fit1$b[5]),exp(fit1$b[6]),
                                2*fit1$b[7],rep(1,ncol(Xmat1)),
                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                                exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                                exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                                2*fit1$b[7],rep(1,ncol(Xmat1)),
                                                                                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                              7+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
        }
        
        
        if(Frailty==0){
          vcov_matrix <- diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                exp(fit1$b[3]),exp(fit1$b[4]),
                                exp(fit1$b[5]),exp(fit1$b[6]),
                                rep(1,ncol(Xmat1)),
                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3)) %*% vcov_matrix %*% diag(c(exp(fit1$b[1]),exp(fit1$b[2]),
                                                                                                exp(fit1$b[3]),exp(fit1$b[4]),
                                                                                                exp(fit1$b[5]),exp(fit1$b[6]),
                                                                                                rep(1,ncol(Xmat1)),
                                                                                                rep(1,ncol(Xmat2)),rep(1,ncol(Xmat3))),
                                                                                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3),
                                                                                              6+ncol(Xmat1)+ncol(Xmat2)+ncol(Xmat3))
        }
        
        
        
        
        terms_object <- terms(formula)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_01 <- ifelse(Frailty == 1, 8, 7)
        
        if(length(covariates_without_cluster01)>0){
          
          covariates_transitions_01 <- fit1$b[debut_reg_01:(debut_reg_01 + ncol(Xmat1)-1)]
          
          cat("Transition 0 -> 1:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat1)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat1))){
            
            cov_name <- colnames(Xmat1)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_01[i], 
                        exp(covariates_transitions_01[i]),
                        sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), 
                        covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]),
                        ifelse(covariates_transitions_01[i] > 0,
                               2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_01[i] / sqrt(vcov_matrix[(debut_reg_01 + i - 1),(debut_reg_01 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          if (length(global_chisq.01) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.01)[dof_chisq.01 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.01[factor_name],
                            dof_chisq.01[factor_name],
                            p.global_chisq.01[factor_name]
                ))
              }
              cat("\n") 
            }
          }
          
        }
        
        
        
        
        
        terms_object <- terms(formula.terminalEvent)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_02 <- debut_reg_01 + ncol(Xmat1)
        if(length(covariates_without_cluster02)>0){
          
          covariates_transitions_02 <- fit1$b[debut_reg_02:(debut_reg_02 +ncol(Xmat2)-1)]
          
          cat("Transition 0 -> 2:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat2)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat2))){
            
            cov_name <- colnames(Xmat2)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_02[i], 
                        exp(covariates_transitions_02[i]),
                        sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), 
                        covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]),
                        ifelse(covariates_transitions_02[i] > 0,
                               2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_02[i] / sqrt(vcov_matrix[(debut_reg_02 + i - 1),(debut_reg_02 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          if (length(global_chisq.02) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.02)[dof_chisq.02 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.02[factor_name],
                            dof_chisq.02[factor_name],
                            p.global_chisq.02[factor_name]
                ))
              }
              cat("\n") 
            }
          }
          
        }
        
        
        
        
        terms_object <- terms(formula.terminalEvent)
        covariates <- attr(terms_object, "term.labels")
        covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
        debut_reg_12 <- debut_reg_02 + ncol(Xmat2)
        if(length(covariates_without_cluster12)>0){
          
          covariates_transitions_12 <- fit1$b[debut_reg_12:(debut_reg_12 + ncol(Xmat3)-1)]
          
          cat("Transition 1 -> 2:\n")
          cat("------------\n")
          
          max_cov_name_length <- max(nchar(colnames(Xmat3)))
          
          
          cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                      max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
          
         
          for(i in 1:length(colnames(Xmat3))){
            
            cov_name <- colnames(Xmat3)[i]
            
            
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n", 
                        max_cov_name_length, cov_name, 
                        covariates_transitions_12[i], 
                        exp(covariates_transitions_12[i]),
                        sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), 
                        covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]),
                        ifelse(covariates_transitions_12[i] > 0,
                               2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = FALSE)),
                               2 * (pnorm(covariates_transitions_12[i] / sqrt(vcov_matrix[(debut_reg_12 + i - 1),(debut_reg_12 + i - 1)]), mean = 0, sd = 1, lower.tail = TRUE)))))
          }
          
          cat("\n")
          
          if (length(global_chisq.12) > 0) {
            
            factors_with_mult_dof <- names(dof_chisq.12)[dof_chisq.12 > 1]
            
            if (length(factors_with_mult_dof) > 0) { 
              cat(sprintf(" %-*s %12s %12s %12s \n",
                          max_cov_name_length, "", "chisq", "df", "global p"))
              
              for (factor_name in factors_with_mult_dof) {
                cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                            max_cov_name_length, factor_name,
                            global_chisq.12[factor_name],
                            dof_chisq.12[factor_name],
                            p.global_chisq.12[factor_name]
                ))
              }
              cat("\n") 
            }
          }
        }
        
        
        # Frailty Parameters
        if(Frailty == 1){
          cat("Frailty Parameters:\n")
          cat("------------\n")
          cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", ((fit1$b[7])^2), sqrt(vcov_matrix[7,7]), 
                      1- pnorm((((fit1$b[7])^2)) / sqrt(vcov_matrix[7,7]))))
          cat("\n")
        }
        
        # Weibull baseline hazard parameters
        cat("Scales and shapes of the Weibull baseline hazard\n")
        cat("------------\n")
        cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
        
       
        cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[1]), sqrt(vcov_matrix[1,1]),
                    exp(fit1$b[2]), sqrt(vcov_matrix[2,2])))
        cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[3]), sqrt(vcov_matrix[3,3]),
                    exp(fit1$b[4]), sqrt(vcov_matrix[4,4])))
        cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", exp(fit1$b[5]), sqrt(vcov_matrix[5,5]),
                    exp(fit1$b[6]), sqrt(vcov_matrix[6,6])))
        
        cat("\n")
        
        cat("The expression of the Weibull hazard function is: \n")
        cat("           'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' \n" )
        cat("The expression of the Weibull survival function is: \n")
        cat("           'S(t) = exp[- (t/scale)^shape]' \n" )
        
        cat("\n")
        
        # Marginal log-likelihood and AIC
        cat("Marginal log-likelihood = ", fit1$fn.value, "\n")
        cat("AIC = Akaike information Criterion = ", (1/nrow(data)) * (length(fit1$b) - fit1$fn.value), "\n")
        cat("           'AIC = (1/n)[np - l(.)]' \n")
        
        cat("\n")
        
        # Number of subjects and transitions
        cat("Number of subjects = ", nrow(data), "\n")
        cat("             0 -> 1 = ", length(which(delta1 == 1)), "\n")
        cat("             0 -> 2 = ", length(which(delta1 == 0 & delta2 == 1)), "\n")
        cat("             1 -> 2 = ", length(which(delta1 == 1 & delta2 == 1)), "\n")
        cat("Lost to follow-up  = ", length(which(delta1 == 0 & delta2 == 0)), "\n")
        cat("\n")
        
        cat("Number of iterations: ", fit1$ni, "\n")
        cat("Convergence criteria:\n")
        cat("------------\n")
        cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
        cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
        cat("Maximum number of iterations reached, parameters in partialH were dropped from Hessian to define the relative distance to optimum\n")
        if (fit1$ca > LIMparam) {
          cat("  Parameter criterion not reached. Threshold: ", LIMparam, "\n")
        }
        
        if (fit1$cb > LIMlogl) {
          cat("  Function criterion not reached. Threshold: ", LIMlogl, "\n")
        }
        
        if (fit1$rdm > LIMderiv) {
          cat("  Relative distance to optimum criterion not reached. Threshold: ", LIMderiv, "\n")
        }
        
      }
      
    }
  }
  
  
  
  
  
  
  Frailty <- ifelse(Frailty==1,TRUE,FALSE)
  
  
  
  
  
  cat("\n")
  
  
  
  
  if(Frailty==TRUE){
    
    if(nvar>0){
      
      if(length(covfactors)==0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(length(covfactors)>0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,
                      names.factor.01=covfactors1,names.factor.02=covfactors2,names.factor.12=covfactors3,
                      global_chisq.01=global_chisq.01,global_chisq.02=global_chisq.02,
                      global_chisq.12=global_chisq.12,
                      p.global_chisq.01=p.global_chisq.01,p.global_chisq.02=p.global_chisq.02,
                      p.global_chisq.12=p.global_chisq.12,
                      dof_chisq.01=dof_chisq.01,dof_chisq.02=dof_chisq.02,dof_chisq.12=dof_chisq.12,
                      global_chisq.test.01=global_chisq.test.01,global_chisq.test.02=global_chisq.test.02,
                      global_chisq.test.12=global_chisq.test.12,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    if(nvar==0){
      
      if(length(covfactors)==0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(length(covfactors)>0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,
                      names.factor.01=covfactors1,names.factor.02=covfactors2,names.factor.12=covfactors3,
                      global_chisq.01=global_chisq.01,global_chisq.02=global_chisq.02,
                      global_chisq.12=global_chisq.12,
                      p.global_chisq.01=p.global_chisq.01,p.global_chisq.02=p.global_chisq.02,
                      p.global_chisq.12=p.global_chisq.12,
                      dof_chisq.01=dof_chisq.01,dof_chisq.02=dof_chisq.02,dof_chisq.12=dof_chisq.12,
                      global_chisq.test.01=global_chisq.test.01,global_chisq.test.02=global_chisq.test.02,
                      global_chisq.test.12=global_chisq.test.12,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    
    
  }
  
  if(Frailty==FALSE){
    
    if(nvar>0){
      if(length(covfactors)==0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv
        )
      }
      
      if(length(covfactors)>0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,
                      names.factor.01=covfactors1,names.factor.02=covfactors2,names.factor.12=covfactors3,
                      global_chisq.01=global_chisq.01,global_chisq.02=global_chisq.02,
                      global_chisq.12=global_chisq.12,
                      p.global_chisq.01=p.global_chisq.01,p.global_chisq.02=p.global_chisq.02,
                      p.global_chisq.12=p.global_chisq.12,
                      dof_chisq.01=dof_chisq.01,dof_chisq.02=dof_chisq.02,dof_chisq.12=dof_chisq.12,
                      global_chisq.test.01=global_chisq.test.01,global_chisq.test.02=global_chisq.test.02,
                      global_chisq.test.12=global_chisq.test.12,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    
    if(nvar==0){
      
      
      if(length(covfactors)==0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      if(length(covfactors)>0){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x01=x01,x02=x02,x12=x12,
                      lam01=lam01,lam02=lam02,lam12=lam12,
                      surv01=surv01, surv02=surv02, surv12=surv12,
                      median.01=median.01, median.02=median.02, median.12=median.12,
                      linear.pred01=linear.pred01, linear.pred02=linear.pred02,
                      linear.pred12=linear.pred12,partialH=partialH,
                      names.factor.01=covfactors1,names.factor.02=covfactors2,names.factor.12=covfactors3,
                      global_chisq.01=global_chisq.01,global_chisq.02=global_chisq.02,
                      global_chisq.12=global_chisq.12,
                      p.global_chisq.01=p.global_chisq.01,p.global_chisq.02=p.global_chisq.02,
                      p.global_chisq.12=p.global_chisq.12,
                      dof_chisq.01=dof_chisq.01,dof_chisq.02=dof_chisq.02,dof_chisq.12=dof_chisq.12,
                      global_chisq.test.01=global_chisq.test.01,global_chisq.test.02=global_chisq.test.02,
                      global_chisq.test.12=global_chisq.test.12,data=data_init,model=model
                      ,formula=formula,formula.terminalEvent=formula.terminalEvent
                      ,Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                      delta1=delta1,delta2=delta2,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    
  }
  
  
  class(value)= "frailtyIllnessDeath"
  
  return(invisible(value))
}







