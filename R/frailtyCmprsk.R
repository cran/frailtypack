#' Fit a Weibull Competing Risks Model with Optional Shared Frailty
#'
#' @description Fit a Weibull competing risks model with shared gamma frailty between
#' all transitions (0->1,0->2,...,0->k). Handles left-truncated and right-censored data.
#' The model considers transitions from an initial state (0) to (k) competing
#' absorbing states.
#'
#' \if{html}{
#' \figure{CMPRSKSCHEME.png}{options: width="329"}
#' }
#' 
#'\if{latex}{
#' \figure{CMPRSKSCHEME.png}{options: width=8.7cm}
#' }
#'
#' @details Let \eqn{T} be the time to event and \eqn{L \in \{1,2,...,k\}}  the indicator of
#' the cause of the event.
#'
#' The cause-specific hazard rate for cause \eqn{l \in \{1,2,...,k\}} is:
#' \deqn{
#'   \lambda_{l}(t) = \lim_{\Delta t \to 0^+} \frac{\mathbb{P}(t \leq T \leq t + \Delta t, L=l \mid T \geq t)}{\Delta t}
#' }

#'
#' A proportional hazards model with a shared frailty term \eqn{\omega_i} is assumed
#' for each transition within group \eqn{i}. For the \eqn{j^{th}} subject
#' (\eqn{j=1,...,n_i}) in the \eqn{i^{th}} group (\eqn{i=1,...,G}), the \eqn{l^{th}}
#' (\eqn{l=1,...,k})  transition intensity is defined as follows:
#' 
#' \deqn{
#'   \lambda_{l}^{ij}(t |\omega_i,X_{l}^{ij}) = \lambda_{0l}(t) \omega_i \exp(\beta_l^{T} X_{l}^{ij})
#' }
#' where \eqn{\omega_i \sim\Gamma(\frac{1}{\theta},\frac{1}{\theta})} with
#' \eqn{\bold{E}(\omega_i)=1} and \eqn{\bold{Var}(\omega_i)=\theta}.
#'
#' \eqn{\omega_i} is the frailty term for the \eqn{i^{th}} group.
#' For subject-specific frailties, use \code{cluster(id)} where id is unique (\eqn{n_i=1}).


#' \eqn{\beta_l} (\eqn{l=1,...,k}) is the vector 
#' of time fixed regression coefficients for the transition 0->l.

#' \eqn{X_{l}^{ij}} (\eqn{l=1,...,k}) is the vector 
#' of time fixed covariates for the \eqn{j^{th}} subject in the \eqn{i^{th}} group for the 
#' transition 0->l.

#' \eqn{\lambda_{0l}(.)} (\eqn{l=1,...,k}) is the 
#' baseline hazard function for the transition 0->l.

#'

#' The Weibull baseline hazard parameterization is:
#' \deqn{\lambda(t) = \frac{\gamma}{\lambda^\gamma} \cdot t^{\gamma - 1}} 


#' where \eqn{\lambda} is the scale parameter and \eqn{\gamma} the shape parameter
#'
#' 
#' @usage frailtyCmprsk(formulas, data, maxit = 300,
#' init.B, init.Theta, init.hazard.weib,
#' LIMparam = 1e-3, LIMlogl = 1e-3, LIMderiv = 1e-3,
#' x0, print.info = FALSE, print.result = TRUE,
#' partialH, blinding = TRUE)
#'
#' @param formulas A list of formula objects. The first formula must include a 
#'   response on the left-hand side of a \code{~} operator. The response must be a
#'   survival object as returned by the \code{Surv} function (e.g., \code{survival::Surv}).
#'   The argument \code{type = "mstate"} must be specified within the
#'   \code{Surv} function. The status indicator in the \code{Surv} object should be:
#'   0 for right-censoring, 1 for the first competing event, 2 for the second, ..., and \code{k} for the
#'   \code{k}th competing event.
#'   Covariates for the transition from 0 to 1 are specified on the right-hand side (RHS) of the first formula.
#'   The remaining elements of the list should be one-sided formulas (e.g., \code{~ var1 + var2}),
#'   used only to specify the covariates. The second formula corresponds to the transition 0 -> 2,
#'   the third to 0 -> 3, and so on, up to the \code{k}th formula for transition 0 -> \code{k}.
#'   Left-truncation is supported and should be specified using the three-argument
#'   \code{Surv(time1, time2, status)} notation.
#'   Shared frailty can be specified via \code{cluster(group_variable)} on the RHS of the first formula only.
#'   It should not be included in the other formulas.

#' @param data A \code{data.frame} containing the variables named in \code{formulas}.
#' @param maxit Maximum number of iterations for the Marquardt optimization algorithm.
#'   Default is 300.
#' @param init.B Optional. A vector of initial values for the regression coefficients.
#'   It should contain the coefficients for each transition from 0->1, ..., up to 0->\code{k}
#'   (i.e., \eqn{\beta_1}, \eqn{\beta_2}, ..., \eqn{\beta_k}), concatenated in order.
#'   The total length must match the total number of covariates specified across all formulas.
#'   If omitted, the regression coefficients are initialized from default values. 
#'   These defaults are obtained by fitting \code{k} independent Weibull proportional 
#'   hazards models (without frailty), one for each transition.
#' @param init.Theta Optional. Initial value for the frailty variance \eqn{\theta}.
#'   Default is 0.1. This parameter is only used if \code{cluster()} is present
#'   in \code{formulas}.
#' @param init.hazard.weib Optional. A vector of initial values for the Weibull baseline
#'   hazard parameters. It must be of length \code{2 * k}, with the values ordered as:
#'   scale(0->1), shape(0->1), scale(0->2), shape(0->2), ..., scale(0->k), shape(0->k).
#'   If omitted, the baseline hazard parameters are initialized from default values. 
#'   These defaults are obtained by fitting \code{k} independent Weibull proportional 
#'   hazards models (without frailty), one for each transition.
#' @param LIMparam Convergence threshold for the parameters based on the maximum
#'   absolute difference between successive iterations (\eqn{10^{-3}} by default).
#' @param LIMlogl Convergence threshold for the log-likelihood based on the absolute
#'   difference between successive iterations (\eqn{10^{-3}} by default).
#' @param LIMderiv Convergence threshold based on the relative distance to the optimum
#'   (related to gradient and Hessian) (\eqn{10^{-3}} by default). See Details.
#' @param x0 Optional. A list of numeric vectors, where each vector specifies the time points 
#'   at which to compute the baseline hazard and survival functions for a given transition. 
#'   The order must follow the transitions: 0->1, 0->2, ..., 0->k. If not provided, defaults to 
#'   a list where each element is a sequence of 99 time points from 0 to the maximum observed 
#'   time for the corresponding transition.
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
#'
#' @return An object of class 'frailtyCmprsk' containing:
#' \describe{
#'   \item{b}{Vector of the estimated parameters. Order is:
#'   (scale(0->1), shape(0->1), scale(0->2), shape(0->2),
#'    ..., scale(0->k), shape(0->k),\eqn{\hat {\theta}} (if frailty), \eqn{\hat{\beta}_1}, \eqn{\hat{\beta}_2},...,\eqn{\hat{\beta}_k}).}
#'   \item{call}{The matched function call.}
#'   \item{coef}{Vector of estimated regression coefficients.}
#'   \item{loglik}{The marginal log-likelihood value at the final parameter estimates.}
#'   \item{grad}{Gradient vector of the log-likelihood at the final parameter estimates.}
#'   \item{n}{The number of observations used in the fit.}
#'   \item{n.events}{Vector containing the number of observed events: count for 0->1, count for 0->2,..., count for 0->k, count for censoring.}
#'   \item{n.iter}{Number of iterations.}
#'   \item{vcov}{Variance-covariance matrix for the parameters listed in \code{b}.}
#'   \item{npar}{Total number of estimated parameters.}
#'   \item{nvar}{Total number of regression coefficients.}
#'   \item{shape.weib}{Vector of estimated Weibull baseline shape parameters (shape(0->1),shape(0->2),...,shape(0->k)).}
#'   \item{scale.weib}{Vector of estimated Weibull baseline scale parameters (scale(0->1),scale(0->2),...,scale(0->k)).}
#'   \item{crit}{Convergence status code: 1=converged, 2=maximum iterations reached, 3=converged using partial Hessian, 4=the algorithm encountered a problem in the loglikelihood computation.}
#'   \item{Frailty}{Logical. \code{TRUE} if a model with shared frailty (\code{cluster(.)}) was fitted.}
#'   \item{beta_p.value}{Vector of p-values from Wald tests for the regression coefficients in \code{coef}.}
#'   \item{AIC}{Akaike Information Criterion, calculated as \eqn{AIC=\frac{1}{n}(np - l(.))}, where np is the number of parameters and l is the log-likelihood.}
#'   \item{x0}{List of numeric vectors, where each vector contains the time points used for calculating the baseline functions for each transition.}
#'   \item{lam0}{List of matrices, where each matrix contains baseline hazard estimates and 95\% confidence intervals at the corresponding time points in \code{x0}, for each transition.}
#'   \item{surv0}{List of matrices, where each matrix contains baseline survival estimates and 95\% confidence intervals at the corresponding time points in \code{x0}, for each transition.}
#'   \item{medians}{List of matrices, where each matrix contains the estimated median baseline survival time and its 95\% confidence interval for each transition.}
#'   \item{linear.pred}{List of numeric vectors, where each vector contains to the linear predictors for transition 0->l (\eqn{l=1,...,k}). 
#'         For non-frailty models, each element is of the form \eqn{\hat{\beta}_l^{t}X_l}. 
#'         For frailty models, it includes the estimated log-frailty: \eqn{\hat{\beta}_l^{t}X_l + \log(\hat{\omega}_i)}.}
#'   \item{names.factor}{List of character vectors, where each vector contains the factor covariates included in the model for the corresponding transition.}
#'   \item{global_chisq}{List of numeric vectors, where each vector contains the chi-squared statistics from global Wald tests for factor variables for the corresponding transition.}
#'   \item{dof_chisq}{List of integer vectors, where each vector contains the degrees of freedom for the global Wald tests for the corresponding transition.}
#'   \item{p.global_chisq}{List of numeric vectors, where each vector contains the p-values for the global Wald tests for the corresponding transition.}
#'   \item{global_chisq.test}{Binary vector of length k indicating whether any global factor tests were performed for each transition.}

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
#' 
#' 
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
#' Marquardt, D. W. (1963). An algorithm for least-squares estimation of nonlinear
#' parameters. \emph{SIAM Journal on Applied Mathematics}, 11(2), 431-441.
#'
#' Liquet, B., Timsit, J. F., & Rondeau, V. (2012). Investigating hospital
#' heterogeneity with a multi-state frailty model: application to nosocomial
#' pneumonia disease in intensive care units. \emph{BMC Medical Research
#' Methodology}, 12(1), 1-14.

#' @examples

#' \donttest{
#' 
#' 
#' data(CPRSKbmtcrr)
#'
#' ##-- Weibull competing risks model with 
#' ##   group frailty shared between transitions --##
#'
#' ModCmprsk_Group <- frailtyCmprsk(
#'   formulas = list(
#'     Surv(observed_time, Status, type = "mstate") ~ cluster(group) + Sex,
#'     ~ Sex
#'   ),
#'   data        = CPRSKbmtcrr,
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' ##-- Weibull competing risks model with subject-specific 
#' ##   frailty shared between transitions --##
#'
#' ModCmprsk_Subject <- frailtyCmprsk(
#'   formulas = list(
#'     Surv(observed_time, Status, type = "mstate") ~ cluster(id) + Sex,
#'     ~ Sex
#'   ),
#'   data        = CPRSKbmtcrr,
#'   print.info  = FALSE,
#'   maxit       = 100
#' )
#'
#' ##--- Simple Weibull competing risks model with left truncation ---##
#'
#' ModCmprsk_LeftTrunc <- frailtyCmprsk(
#'   formulas = list(
#'     Surv(Age, observed_time, Status, type = "mstate") ~ Source,
#'     ~ Sex
#'   ),
#'   data       = CPRSKbmtcrr,
#'   print.info = FALSE
#' )
#'
#' ##--- Simple Weibull competing risks model with a factor and 
#' ##    no covariates for the first competing event (left truncation) ---##
#'
#' ModCmprsk_Factor_LeftTrunc <- frailtyCmprsk(
#'   formulas = list(
#'     Surv(Age, observed_time, Status, type = "mstate") ~ Source,
#'     ~ factor(Phase)
#'   ),
#'   data       = CPRSKbmtcrr,
#'   print.info = FALSE
#' )  }
#'
#' @import marqLevAlg 
#' @importFrom stats delete.response model.response reformulate





















#' @export
"frailtyCmprsk" <- function(formulas, 
                            data , maxit = 300,
                            init.B, init.Theta, init.hazard.weib,
                            LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, 
                            x0, print.info=FALSE, print.result=TRUE, 
                            partialH, blinding=TRUE){
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## Log-likelihood illness-death, shared frailty, left-truncated data
  ## Weibull baseline hazards for n competing risks
  logLike.group.weibull.SCR.SM.LT.FRAILTY <- function(para, y, delta, l, group, data, Xmat_list)
  {
    
    delta <- as.numeric(delta)
    
    
    n_events <- max(delta)
    
    delta_indicators <- lapply(1:n_events, function(k) ifelse(delta == 0, 0, ifelse(delta == k, 1, 0)))
    
    
    kappa <- numeric(n_events)
    alpha <- numeric(n_events)
    
    for (k in 1:n_events) {
      kappa[k] <- exp((para[2*k - 1]))  
      alpha[k] <- exp(para[2*k])           
    }
    
    theta    <- (para[2*n_events + 1])^2 
    thetaInv <- 1 / theta
    
    nP.fixed <- sum(sapply(Xmat_list, ncol))
    
    
    beta_start_idx <- 2*n_events + 2
    
    eta_list <- vector("list", n_events)
    current_beta_idx <- beta_start_idx
    for (k in 1:n_events) {
      current_Xmat <- Xmat_list[[k]]
      if (ncol(current_Xmat) == 0) {
        eta_list[[k]] <- rep(0, nrow(current_Xmat))
      } else {
        num_betas_k <- ncol(current_Xmat)
        eta_list[[k]] <- as.vector(current_Xmat %*% para[current_beta_idx:(current_beta_idx + num_betas_k - 1)])
        current_beta_idx <- current_beta_idx + num_betas_k
      }
    }
    
    log_h_star_y_list <- vector("list", n_events)
    q_y_components <- matrix(0, nrow = length(y), ncol = n_events) 
    q_l_components <- matrix(0, nrow = length(l), ncol = n_events) 
    
    for (k in 1:n_events) {
      log_h_star_y_list[[k]] <- log(alpha[k]) - alpha[k]*log(kappa[k]) + (alpha[k] - 1) * log(y) + eta_list[[k]]
      q_y_components[,k] <- (y/kappa[k])^alpha[k] * exp(eta_list[[k]])
      q_l_components[,k] <- (l/kappa[k])^alpha[k] * exp(eta_list[[k]])
    }
    
    
    q.y <- rowSums(q_y_components)
    q.l <- rowSums(q_l_components)
    
    
    loglike <- 0
    
    for(i in unique(group)){
      
      ni <- which(data$group == i)
      
      mi_total <- 0
      for(k in 1:n_events){
        mi_total <- mi_total + sum(delta_indicators[[k]][ni])
      }
      
      for(j in ni){
        for(k in 1:n_events){
          loglike <- loglike + delta_indicators[[k]][j] * log_h_star_y_list[[k]][j]
        }
      }
      
      loglike <- loglike - thetaInv * log(theta)
      
      if(mi_total > 0){
        
        
        for(k_prime in 0:(mi_total - 1)){ 
          loglike <- loglike + log(thetaInv + k_prime)
        }
      }
      
      
      v_group_sum = 0
      for(j in ni){
        v_group_sum <- v_group_sum + q.y[j]
      }
      loglike <- loglike - (mi_total + thetaInv) * log(v_group_sum + thetaInv)
      
      
      lt_group_sum = 0
      for (j in ni){
        lt_group_sum <- lt_group_sum + q.l[j]
      }
      loglike <- loglike + thetaInv * log(1 + theta * lt_group_sum)
      
    }
    
    return(loglike)
  }
  
  
  ################
  ## Log-likelihood illness-death, NO shared frailty, left-truncated data
  ## Weibull baseline hazards for n competing risks
  logLike.weibull.SCR.SM.LT.SANS.FRAILTY <- function(para, y, delta, l, Xmat_list)
  {
    
    delta <- as.numeric(delta)
    
    
    n_events <- max(delta)
    
    
    delta_indicators <- lapply(1:n_events, function(k) ifelse(delta == k, 1, 0))
    
    censored_indicator <- ifelse(delta == 0, 1, 0) 
    
    
    
    
    kappa <- numeric(n_events)
    alpha <- numeric(n_events)
    
    for (k in 1:n_events) {
      kappa[k] <- exp((para[2*k - 1]))  
      alpha[k] <- exp(para[2*k])          
    }
    
    
    beta_start_idx <- 2*n_events + 1 
    
    
    eta_list <- vector("list", n_events)
    current_beta_idx <- beta_start_idx
    for (k in 1:n_events) {
      current_Xmat <- Xmat_list[[k]]
      if (ncol(current_Xmat) == 0) {
        eta_list[[k]] <- rep(0, nrow(current_Xmat))
      } else {
        num_betas_k <- ncol(current_Xmat)
        eta_list[[k]] <- as.vector(current_Xmat %*% para[current_beta_idx:(current_beta_idx + num_betas_k - 1)])
        current_beta_idx <- current_beta_idx + num_betas_k
      }
    }
    
    
    log_h_star_y_list <- vector("list", n_events)
    q_y_components <- matrix(0, nrow = length(y), ncol = n_events) 
    q_l_components <- matrix(0, nrow = length(l), ncol = n_events) 
    
    for (k in 1:n_events) {
      log_h_star_y_list[[k]] <- log(alpha[k]) - alpha[k]*log(kappa[k]) + (alpha[k] - 1) * log(y) + eta_list[[k]]
      q_y_components[,k] <- (y/kappa[k])^alpha[k] * exp(eta_list[[k]])
      q_l_components[,k] <- (l/kappa[k])^alpha[k] * exp(eta_list[[k]])
    }
    
    
    q.y <- rowSums(q_y_components)
    q.l <- rowSums(q_l_components)
    
    
    loglh <- 0
    
    for (j in 1:length(y)) { 
      if (censored_indicator[j] == 1) {
        
        loglh <- loglh + (- q.y[j] + q.l[j])
      } else {
        
        
        event_type_j <- delta[j] 
        if (event_type_j > 0 && event_type_j <= n_events) { 
          loglh <- loglh + log_h_star_y_list[[event_type_j]][j] - q.y[j] + q.l[j]
        }
      }
    }
    
    return(loglh)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Function to calculate confidence bands for a set of times
  
  hazard_and_confidence_bands_cmprsk <- function(t_vec,scale,shape,varcov) {
    
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
  
  survival_and_confidence_bands_cmprsk <- function(t_vec,scale,shape,varcov) {
    
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
  
  
  
  
  
  
  
  
  
  
  
  
  # Function to calculate hazars for a set of times
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
  
  
  
  ### CHECK IF ARGUMENTS ARE OF CORRECT TYPES
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
  
  
  n_competing_events_from_formulas <-  length(formulas)
  
  
  
  
  
  if(!missing(init.Theta)){
    if (!is.numeric(init.Theta) || length(init.Theta) != 1 || init.Theta < 0) {
      stop("'init.Theta' should be a positive numeric value.")
    }
  }
  
  ### CHECKING FORMULAS AND COVARIATES VALIDITY
  
  
  
  if (!is.list(formulas)) {
    stop("'formulas' should be a list of formula objects.")
  }
  
  
  if (length(formulas) == 0) {
    stop("'formulas' list cannot be empty. It must contain at least two formulas.")
  }
  
  
  for (i in seq_along(formulas)) {
    if (!inherits(formulas[[i]], "formula")) {
      stop(paste0("Element ", i, " in 'formulas' is not a formula object."))
    }
  }
  
  
  tryCatch({
    
    mf_first_formula <- model.frame(formulas[[1]], data = data, na.action = na.pass)
    
    response_obj <- model.response(mf_first_formula)
  }, error = function(e) {
    stop(paste0("Could not extract response from the first formula. Ensure all variables are in 'data'. Error: ", e$message))
  })
  
  
  if (!inherits(response_obj, "Surv")) {
    stop("The left-hand side of the first formula in 'formulas' must evaluate to a Surv object.")
  }
  
  
  if (attr(response_obj, "type") != "mright" && attr(response_obj, "type") != "mcounting") {
    stop("The 'Surv' object in the first formula must specify 'type=\"mstate\"'.")
  }
  
  
  
  if (!ncol(response_obj) %in% c(2, 3)) {
    stop("The 'Surv' object for 'mstate' should have 2 or 3 columns (time/status or time1/time2/status).")
  }
  
  
  
  
  
  Xmat_list <- vector("list", length(formulas))
  
  Xmat_list_names <- vector("list", length(formulas))
  
  covariates_raw_labels <- vector("list", length(formulas))
  
  
  factor_covariates_per_event <- vector("list", length(formulas))
  
  
  
  
  
  all_unique_covariate_names <- character(0)
  
  
  
  
  group_variable <- NULL
  
  group_data <- NULL
  
  
  
  
  for (k in seq_along(formulas)) {
    
    
    current_formula <- formulas[[k]]
    
    
    
    
    
    terms_object <- terms(current_formula) 
    
    
    
    
    if (k == 1) {
      
      
      if (attr(terms_object, "response") == 0) {
        
        
        stop(paste0("The first formula must specify the Surv object on its left-hand side."))
        
      }
      
    } else {
      
      
      if (attr(terms_object, "response") == 1) {
        
        
        stop(paste0("Formula ", k , " (", deparse(formulas[[k]]), ") ", "should only specify covariates on its right-hand side (e.g., ~ var1 + var2). Do not include a response on the left-hand side."))
        
      }
      
    }
    
    
    
    
    
    
    
    
    
    
    
    current_covariates_raw_labels <- attr(terms_object, "term.labels")
    
    
    
    cluster_term_indices <- grep("^cluster\\((.*)\\)$", current_covariates_raw_labels)
    
    
    
    
    if(length(cluster_term_indices)!=0){
      covariates_without_cluster <- current_covariates_raw_labels[-cluster_term_indices]
      covariates_with_factor <- grep("^factor\\(.*\\)$", covariates_without_cluster, value = TRUE)
      covariates_with_factor <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates_with_factor))
      if(length(covariates_with_factor)!=0){
        data[[covariates_with_factor]] <- as.factor(data[[covariates_with_factor]])
      }
      covariates_without_cluster <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster))
    }else{
      covariates_without_cluster <- current_covariates_raw_labels
      covariates_with_factor <- grep("^factor\\(.*\\)$", covariates_without_cluster, value = TRUE)
      covariates_with_factor <- unique(gsub("^factor\\((.*)\\)$", "\\1", covariates_with_factor))
      if(length(covariates_with_factor)!=0){
        data[[covariates_with_factor]] <- as.factor(data[[covariates_with_factor]])
      }
      covariates_without_cluster<-unique(gsub("^factor\\((.*)\\)$", "\\1", covariates_without_cluster))
      
      
      
    }
    
    if (length(cluster_term_indices) > 0) { 
      
      if (k == 1) {
        
        if (length(cluster_term_indices) > 1) {
          
          
          stop("Only one 'cluster()' term is allowed in the first formula for a shared frailty model.")
          
        }
        
        
        
        cluster_term_string <- current_covariates_raw_labels[cluster_term_indices[1]]
        
        group_variable_expr <- gsub("^cluster\\((.*)\\)$", "\\1", cluster_term_string)
        
        group_variable <- as.character(group_variable_expr)
        
        
        
        if (!group_variable %in% names(data)) {
          
          
          stop(paste0("The cluster variable '", group_variable, "' specified in formulas[[1]] was not found in the 'data' data frame."))
          
        }
        
        group_data <- data[[group_variable]]
        
        
        
        
        
        
        
        
        
        
        
      } else {
        
        
        stop(paste0("For a SHARED frailty model, the 'cluster()' term must be specified ONLY in the first formula" ))
        
      }
      
    } 
    
    if(k==1){
      if(length(covariates_without_cluster)>0){
        current_formula_for_Xmat <- reformulate(
          
          covariates_without_cluster,
          
          response = current_formula[[2]] 
          
        )
        
        
        
        current_formula_for_Xmat_names <- reformulate(
          
          all.vars(delete.response(terms(current_formula_for_Xmat))),
          
          response = current_formula[[2]] 
          
        )
      }else {
        
        
        current_formula_for_Xmat <- reformulate(
          
          "1",
          
          response = current_formula[[2]] 
          
        )
        
        
        
        current_formula_for_Xmat_names <- all.vars(delete.response(terms(current_formula_for_Xmat)))
        
        
        
        
        
      }
    }else{
      if(length(covariates_without_cluster)>0){
        current_formula_for_Xmat <- reformulate(
          
          covariates_without_cluster,
          
          response = NULL 
          
        )
        
        
        
        current_formula_for_Xmat_names <- reformulate(
          
          all.vars(delete.response(terms(current_formula_for_Xmat))),
          
          response = NULL 
          
        )
      }else {
        
        
        current_formula_for_Xmat <- reformulate(
          
          "1",
          
          response = NULL 
          
        )
        
        
        
        current_formula_for_Xmat_names <- all.vars(delete.response(terms(current_formula_for_Xmat)))
        
        
        
        
        
      }
    }
    
    
    rhs_vars <- all.vars(delete.response(terms(current_formula_for_Xmat)))
    
    
    
    all_unique_covariate_names <- unique(c(all_unique_covariate_names, rhs_vars))
    
    
    if(!is.null(rhs_vars)){
      covariates_raw_labels[[k]] <- rhs_vars
    }
    
    
    
    
    
    
    
    current_event_factor_names <- c()
    
    
    
    
    
    for (var_name in rhs_vars) {
      
      
      
      tryCatch({
        
        var_data <- data[[var_name]] 
        
        if (is.factor(var_data) || is.character(var_data)) {
          
          current_event_factor_names <- c(current_event_factor_names, var_name)
          
        }
        
      }, error = function(e) {
        
        
        warning(paste0("Could not determine class of variable '", var_name, "' for formula ", k, ": ", e$message))
        
      })
      
    }
    if(!is.null(unique(current_event_factor_names))){
      factor_covariates_per_event[[k]] <- unique(current_event_factor_names) 
    }
    
    
    
    
    
    
    if(length(current_formula_for_Xmat_names)!=0){
      tryCatch({
        
        Xmat_list_names[[k]] <- model.matrix(delete.response(terms(current_formula_for_Xmat_names)), data = data)
        
        
        
        if ("(Intercept)" %in% colnames(Xmat_list_names[[k]])) {
          
          Xmat_list_names[[k]] <- Xmat_list_names[[k]][, -which(colnames(Xmat_list_names[[k]]) == "(Intercept)"), drop = FALSE]
          
        }
        
      }, error = function(e) {
        
        
        stop(paste0("Could not create design matrix for formula ", k, ". Ensure all covariates are accessible in 'data'. Error: ", e$message))
        
      })
    }else{
      Xmat_list_names[[k]] <- current_formula_for_Xmat_names
    }
    
    
    
    
    
    tryCatch({
      
      Xmat_list[[k]] <- model.matrix(delete.response(terms(current_formula_for_Xmat)), data = data)
      
      
      
      if ("(Intercept)" %in% colnames(Xmat_list[[k]])) {
        
        Xmat_list[[k]] <- Xmat_list[[k]][, -which(colnames(Xmat_list[[k]]) == "(Intercept)"), drop = FALSE]
        
      }
      
    }, error = function(e) {
      
      
      stop(paste0("Could not create design matrix for formula ", k, ". Ensure all covariates are accessible in 'data'. Error: ", e$message))
      
    })
    
  }
  
  
  
  
  Frailty <- 0
  
  if (!is.null(group_variable)) {
    Frailty <- 1
  }
  
  
  
  
  
  
  
  
  data <- drop_na(data)
  
  
  
  if (!is.null(group_variable)) {
    cluster_variable <- sub("cluster\\((.*?)\\)", "\\1", deparse(group_variable))
    group <- eval(parse(text = paste0("data$", cluster_variable)))
    data$group <- group
  }
  
  if(!missing(init.Theta)){
    
    if(is.null(group_variable)){ 
      stop("init.Theta does not exist in a simple Weibull competing risks model")
    }
    
  }
  
  
  
  
  
  
  lhs <- formulas[[1]][[2]]
  
  if (inherits(lhs, "call") && lhs[[1]] == as.symbol("Surv")) {
    
    surv_terms1 <- as.character(lhs[-1])  
    
    surv_terms1 <- surv_terms1[!grepl("mstate|type\\s*=", surv_terms1)]
    
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
    
    
    formula11 <- paste("Surv(", paste(surv_terms1, collapse = ", "), ", type='mstate')")
    
    formula11 <- eval(parse(text = formula11))
    
  }
  
  
  
  if(!inherits(lhs, "call")) {
    if(inherits(eval(lhs), "Surv")){
      formula11 <- eval(lhs)  
      surv_trunc1 <- attr(eval(lhs), "dimnames")[[2]]
      surv_terms1 <- surv_trunc1
    }
    
  }
  
  
  if(length(surv_trunc1)==3){
    y     <- as.vector(formula11[,"stop"])
    delta <- (as.vector(formula11[,"status"]))
    l <-  as.vector(formula11[,"start"])
  }
  
  if(length(surv_trunc1)==2){
    y     <- as.vector(formula11[,"time"])
    delta <- as.vector(formula11[,"status"])
    l <-  as.vector(rep(0,nrow(data)))
  }
  
  
  n_competing_events_from_formulas <- length(formulas)
  
  valid_status_values <- c(0, 1:n_competing_events_from_formulas)
  
  if (!identical(sort(unique(delta)), valid_status_values)) {
    stop(paste0(
      "The 'status' variable in the Surv object does not exactly match the expected structure based on the number of formulas. ",
      "For a model with ", n_competing_events_from_formulas, " competing events (as implied by your formulas), ",
      "the 'status' variable must contain **only** the unique values: ",
      "0 (for censored), and integers from 1 to ", n_competing_events_from_formulas,
      " (e.g., 1 for the first competing event, 2 for the second, etc.). ",
      "All these expected values must be present, and no other values should appear. ",
      "Detected unique status values: ", paste(sort(unique(delta)), collapse=", "), "."
    ))
  }
  
  
  
  
  delta_indicators <- lapply(1:n_competing_events_from_formulas, function(k) ifelse(delta == k, 1, 0))
  
  
  
  
  
  if (missing(x0)) {
    x0 <- vector("list", n_competing_events_from_formulas)
    for (k in 1:n_competing_events_from_formulas) {
      
      max_time_for_event_k <- max(y[delta_indicators[[k]] == 1])
      x0[[k]] <- round(seq(0, max_time_for_event_k, length.out = 99), 3)
    }
  } 
  
  
  if (!missing(x0)) {
    if (!is.list(x0)) {
      stop("'x0' should be a list of numeric vectors, one vector for each competing event.")
    }
    if (length(x0) != n_competing_events_from_formulas) {
      stop(paste0("Wrong number of elements in 'x0'. Expected ",
                  n_competing_events_from_formulas, " vectors (one for each competing event formula), but received ",
                  length(x0), "."))
    }
    for (k in 1:n_competing_events_from_formulas) {
      current_x0 <- x0[[k]]
      if (!is.numeric(current_x0) || !all(current_x0 >= 0)) {
        stop(paste0("Element ", k, " of 'x0' (for competing event ", k, ") should be a numeric vector with non-negative values."))
      }
      if (length(current_x0) == 0) {
        stop(paste0("Element ", k, " of 'x0' (for competing event ", k, ") cannot be an empty vector."))
      }
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  data_init <- data.frame(y=y,delta=delta,
                          
                          l=l)
  
  
  
  
  
  
  all_beta_names <- c()
  for (k in 1:length(formulas)) {
    current_Xmat_colnames <- colnames(Xmat_list[[k]])
    if (length(current_Xmat_colnames) > 0) {
      
      all_beta_names <- c(all_beta_names, paste0(current_Xmat_colnames, ".event", k))
    }
  }
  
  
  
  
  
  
  
  
  frailty_init_results <- vector("list", n_competing_events_from_formulas)
  initial_scale_weib <- numeric(n_competing_events_from_formulas)
  initial_shape_weib <- numeric(n_competing_events_from_formulas)
  initial_coefs_list <- vector("list", n_competing_events_from_formulas)
  covariates_raw_labels_transormed <- covariates_raw_labels
  
  for (k in 1:n_competing_events_from_formulas) {
    current_Xmat <- Xmat_list[[k]]
    current_delta_indicator <- delta_indicators[[k]]
    
    
    
    
    
    
    temp_data_for_frailtyPenal <- data.frame(
      .l = l,
      .y = y,
      .delta_k = current_delta_indicator
    )
    if (ncol(current_Xmat) > 0) {
      temp_data_for_frailtyPenal <- cbind(temp_data_for_frailtyPenal, data[covariates_raw_labels[[k]]])
      
      is_char_or_factor <- sapply(data[covariates_raw_labels[[k]]], function(x) is.character(x) || is.factor(x))
      
      names_char_or_factor <- names(data[covariates_raw_labels[[k]]])[is_char_or_factor]
      
      indexes_char_or_factor <- which(is_char_or_factor)  
      
      covariates_raw_labels_transormed[[k]][indexes_char_or_factor] <- 
        paste0("factor(", covariates_raw_labels[[k]][indexes_char_or_factor], ")")
      
      covariate_terms_k <- paste(covariates_raw_labels_transormed[[k]], collapse = " + ")
      dynamic_formula_k <- as.formula(paste("Surv(.l, .y, .delta_k) ~", covariate_terms_k))
    } else {
      dynamic_formula_k <- Surv(.l, .y, .delta_k) ~ 1
    }
    
    suppressWarnings({
      frailty_init_results[[k]] <- 
        frailtyPenal(dynamic_formula_k, data = temp_data_for_frailtyPenal,
                     hazard = "Weibull", maxit = maxit, print.times = FALSE)
      
    })
    
    if (frailty_init_results[[k]]$istop == 1) {
      initial_scale_weib[k] <- frailty_init_results[[k]]$scale.weib[1]
      initial_shape_weib[k] <- frailty_init_results[[k]]$shape.weib[1]
      if (ncol(current_Xmat) > 0) {
        initial_coefs_list[[k]] <- frailty_init_results[[k]]$coef
      } else {
        initial_coefs_list[[k]] <- numeric(0) 
      }
    } else {
      
      initial_scale_weib[k] <- 0.1 
      initial_shape_weib[k] <- 0.1 
      if (ncol(current_Xmat) > 0) {
        initial_coefs_list[[k]] <- rep(0.1, ncol(current_Xmat))
        names(initial_coefs_list[[k]]) <- colnames(Xmat_list_names[[k]])
      }
      else {
        initial_coefs_list[[k]] <- numeric(0) 
        names(initial_coefs_list[[k]]) <- colnames(Xmat_list_names[[k]])
      }
    }
    
    
    
    
    weibull_start_vals <- c()
    for (k in 1:n_competing_events_from_formulas) {
      weibull_start_vals <- c(weibull_start_vals,
                              log(initial_scale_weib[k]), 
                              log(initial_shape_weib[k]))      
    }
    
    
  }
  
  theta_start_val <- NULL
  if (Frailty == 1) {
    
    theta_start_val <- sqrt(0.1) 
  }
  
  
  beta_start_vals <- unlist(initial_coefs_list)
  
  if(Frailty==1){
    startVals <- c(weibull_start_vals, theta_start_val, beta_start_vals)
    
  }else{
    startVals <- c(weibull_start_vals, beta_start_vals)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if (!missing(init.hazard.weib)) {
    
    expected_hazard_params_count <- n_competing_events_from_formulas * 2
    
    if (!is.atomic(init.hazard.weib) || length(init.hazard.weib) != expected_hazard_params_count) {
      stop(paste0(
        "'init.hazard.weib' must be a vector of size 2 * ", n_competing_events_from_formulas,
        " ( a scale and a shape for each of the ", n_competing_events_from_formulas, " events), ",
        "with the order: scale1, shape1, scale2, shape2, ... . "
      ))
    }
    
    
    if (any(sapply(init.hazard.weib, function(x) !(is.numeric(x) && x > 0) ))) {
      stop("Each element of 'init.hazard.weib' must be a positive real number.")
    }
    
    for (i in 1:length(init.hazard.weib)) {
    
        
          startVals[i] <- log(init.hazard.weib[i])
       
    }
  }
  
  
  
  
  if (!missing(init.B)) {
    
    if(any(sapply(init.B, function(x) !is.numeric(x)))) {
      stop("init.B should be a vector of real numbers.")
    }
    
    total_expected_betas <- sum(sapply(Xmat_list, ncol))
    
    if (!is.atomic(init.B) || length(init.B) != total_expected_betas) {
      stop(paste0(
        "'init.B' must be a vector of size ", total_expected_betas,
        " (representing the regression coefficients for all covariates across all transitions). ",
        "The order of elements should be: coefficients for event 1, then coefficients for event 2, and so on. ",
        "For each event, coefficients must follow the order of covariates as they appear in the corresponding formula."
      ))
    }
    
    
    
    
    beta_start_index_in_startVals <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1 
    
    current_beta_idx_in_initB <- 1 
    current_beta_idx_in_startVals <- beta_start_index_in_startVals
    
    for (k in 1:n_competing_events_from_formulas) {
      num_betas_k <- ncol(Xmat_list[[k]])
      if (num_betas_k > 0) {
        for (i in 1:num_betas_k) {
          init_B_val <- init.B[current_beta_idx_in_initB]
          if (is.numeric(init_B_val)) {
            startVals[current_beta_idx_in_startVals] <- init_B_val
          } 
          current_beta_idx_in_initB <- current_beta_idx_in_initB + 1
          current_beta_idx_in_startVals <- current_beta_idx_in_startVals + 1
        }
      }
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
  
  
  trunc <- ifelse(length(surv_terms1)==3,TRUE,FALSE)
  
  
  
  
  if(Frailty==1){
    
    
    cat("Be patient the program is computing...\n")
    logLike <- function(p) logLike.group.weibull.SCR.SM.LT.FRAILTY(p, y=y, delta=delta , l=l,
                                                                   data=data,group=group,Xmat_list=Xmat_list)
    error <- NULL
    tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                 ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                 epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
             error = function(e) {
               error <<- e$message
             })
    
    
  }
  
  
  
  
  if(Frailty==0){
    
    
    cat("Be patient the program is computing...\n")
    logLike <- function(p) logLike.weibull.SCR.SM.LT.SANS.FRAILTY(p, y=y, delta=delta , l=l,
                                                                  Xmat_list=Xmat_list)
    error <- NULL
    tryCatch({fit1 <- marqLevAlg(b=startVals,fn=logLike,minimize=FALSE
                                 ,blinding=blinding, maxiter = maxit,print.info = print.info,
                                 epsa = LIMparam,epsb = LIMlogl,epsd=LIMderiv,partialH = partialH)},
             error = function(e) {
               error <<- e$message
             })
    
    
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
  
  
  
  if(Frailty==1){
    
    cat("Competing risks model with shared frailty between transitions\n")
    cat("Using Weibull baseline hazard functions \n")
    if(length(surv_terms1)==3 ){
      cat("Left truncation structure is used\n")
    }
    
    
    cat("\n")
  }
  
  if(Frailty==0){
    
    cat("Competing risks model\n")
    cat("Using Weibull baseline hazard functions \n")
    if(length(surv_terms1)==3 ){
      cat("Left truncation structure is used\n")
    }
    
    
    cat("\n")
    
  }
  
  # DELTA METHOD FOR VCOV MATRIX 
  
  v <- fit1$v
  p <- length(fit1$b)
  vcov_matrix <- matrix(0, p, p)
  vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
  
  vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
  
  
  jacobian_diag_elements <- c()
  current_b_idx <- 1 
  
  for (k in 1:n_competing_events_from_formulas) {
    jacobian_diag_elements <- c(jacobian_diag_elements, exp(fit1$b[current_b_idx]))
    current_b_idx <- current_b_idx + 1
    
    jacobian_diag_elements <- c(jacobian_diag_elements,  exp(fit1$b[current_b_idx]))
    current_b_idx <- current_b_idx + 1
  }
  
  if (Frailty == 1) {
    jacobian_diag_elements <- c(jacobian_diag_elements, 2 * fit1$b[current_b_idx])
    current_b_idx <- current_b_idx + 1
  }
  
  total_beta_count <- sum(sapply(Xmat_list, ncol))
  if (total_beta_count > 0) {
    jacobian_diag_elements <- c(jacobian_diag_elements, rep(1, total_beta_count))
  }
  
  
  ## DELTA METHOD
  jacobian_matrix <- diag(jacobian_diag_elements)
  vcov_matrix <- jacobian_matrix %*% vcov_matrix %*% jacobian_matrix
  
  
  
  
  
  
  
  # --- Global Chi-Squared Tests for Factors ---
  global_chisq <- vector("list", n_competing_events_from_formulas)
  dof_chisq <- vector("list", n_competing_events_from_formulas)
  p.global_chisq <- vector("list", n_competing_events_from_formulas)
  global_chisq.test <- numeric(n_competing_events_from_formulas)
  names.factor <- vector("list", n_competing_events_from_formulas)
  
  
  beta_start_index_in_b <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1
  
  current_beta_idx_in_b <- beta_start_index_in_b
  current_vcov_idx_for_factors <- beta_start_index_in_b
  
  
  for (k in 1:n_competing_events_from_formulas) {
    current_covfactors <- factor_covariates_per_event[[k]]
    
    current_Xmat <- Xmat_list[[k]]
    
    num_betas_in_current_event <- ncol(current_Xmat)
    
    if (length(current_covfactors) > 0) {
      global_chisq_k <- c()
      dof_chisq_k <- c()
      p.global_chisq_k <- c()
      
      if (num_betas_in_current_event > 0) {
        vcov_coef_k <- vcov_matrix[current_vcov_idx_for_factors:(current_vcov_idx_for_factors + num_betas_in_current_event - 1),
                                   current_vcov_idx_for_factors:(current_vcov_idx_for_factors + num_betas_in_current_event - 1), drop = FALSE]
        
        coef_k <- fit1$b[current_beta_idx_in_b:(current_beta_idx_in_b + num_betas_in_current_event - 1)]
        
        
        for (factor_name in current_covfactors) {
          clean_factor_name <- gsub("^factor\\((.*)\\)$", "\\1", factor_name)
          
          
          cols_for_factor_in_Xmat_logical <- startsWith(colnames(current_Xmat), clean_factor_name)
          cols_for_factor_in_Xmat <- which(cols_for_factor_in_Xmat_logical)
          
          dof <- length(cols_for_factor_in_Xmat)
          
          if (dof > 0) {
            factor_coefs <- coef_k[cols_for_factor_in_Xmat]
            factor_vcov <- vcov_coef_k[cols_for_factor_in_Xmat, cols_for_factor_in_Xmat, drop = FALSE]
            
            W <- NA 
            p <- NA 
            
            if (dof == 1) { 
              W <- (factor_coefs^2) / factor_vcov[1,1]
            } else { 
              if (rcond(factor_vcov) < .Machine$double.eps * nrow(factor_vcov)) { 
                warning(paste0("Singular variance-covariance matrix for factor '", factor_name, "' in event ", k, ". Skipping global test for this factor."))
              } else {
                W <- t(factor_coefs) %*% solve(factor_vcov) %*% factor_coefs
              }
            }
            
            if (!is.na(W)) {
              p <- 1 - pchisq(W, df = dof)
            }
            
            global_chisq_k <- c(global_chisq_k, W)
            dof_chisq_k <- c(dof_chisq_k, dof)
            p.global_chisq_k <- c(p.global_chisq_k, p)
          }
        } 
      } 
      
      if (length(global_chisq_k) > 0) {
        names.factor[[k]] <- current_covfactors
        names(global_chisq_k) <- current_covfactors
        names(dof_chisq_k) <- current_covfactors
        names(p.global_chisq_k) <- current_covfactors
        global_chisq[[k]] <- global_chisq_k
        dof_chisq[[k]] <- dof_chisq_k
        p.global_chisq[[k]] <- p.global_chisq_k
        global_chisq.test[k] <- 1 
      } else {
        
        global_chisq.test[k] <- 0
      }
    } else { 
      
      global_chisq.test[k] <- 0
    }
    
    
    current_beta_idx_in_b <- current_beta_idx_in_b + num_betas_in_current_event
    current_vcov_idx_for_factors <- current_vcov_idx_for_factors + num_betas_in_current_event
  } 
  
  names(global_chisq) <- paste0("event", 1:n_competing_events_from_formulas)
  names(dof_chisq) <- paste0("event", 1:n_competing_events_from_formulas)
  names(p.global_chisq) <- paste0("event", 1:n_competing_events_from_formulas)
  names(global_chisq.test) <- paste0("event", 1:n_competing_events_from_formulas)
  names(names.factor) <- paste0("event", 1:n_competing_events_from_formulas)
  
  
  
  
  
  
  
  
  if(direct==TRUE){
    
    
    
    
    
    
    
    
    current_b_idx_for_betas <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1 
    current_vcov_idx_for_betas <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1 
    
    for (k in 1:n_competing_events_from_formulas) {
      current_Xmat <- Xmat_list[[k]]
      num_betas_k <- ncol(current_Xmat)
      
      if (num_betas_k > 0) {
        cat(paste0("Transition 0 -> ", k, ":\n"))
        cat("------------\n")
        
        betas_k <- fit1$b[current_b_idx_for_betas:(current_b_idx_for_betas + num_betas_k - 1)]
        se_betas_k <- sqrt(diag(vcov_matrix)[current_vcov_idx_for_betas:(current_vcov_idx_for_betas + num_betas_k - 1)])
        
        z_values_k <- betas_k / se_betas_k
        p_values_k <- 2 * pnorm(abs(z_values_k), lower.tail = FALSE) # Two-sided p-value
        
        max_cov_name_length <- max(nchar(colnames(current_Xmat)), 10) # Min 10 for header alignment
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n",
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        for (i in 1:num_betas_k) {
          cov_name <- colnames(current_Xmat)[i]
          cat(sprintf(" %-*s %12.5f %12.5f %12.5f %12.5f %12.5f\n",
                      max_cov_name_length, cov_name,
                      betas_k[i],
                      exp(betas_k[i]),
                      se_betas_k[i],
                      z_values_k[i],
                      p_values_k[i]))
        }
        cat("\n")
        
        # GLOBAL TESTS PRINT
        if (length(global_chisq[[k]]) > 0) {
          factors_with_mult_dof <- names(dof_chisq[[k]])[dof_chisq[[k]] > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          global_chisq[[k]][factor_name],
                          dof_chisq[[k]][factor_name],
                          p.global_chisq[[k]][factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
        current_b_idx_for_betas <- current_b_idx_for_betas + num_betas_k
        current_vcov_idx_for_betas <- current_vcov_idx_for_betas + num_betas_k
      }
    }
    
    
    if (Frailty == 1) {
      cat("Frailty Parameters:\n")
      cat("------------\n")
      theta_index_in_b <- length(weibull_start_vals) + 1 
      
      theta_val <- (fit1$b[theta_index_in_b])^2
      
      
      se_theta <- sqrt(vcov_matrix[theta_index_in_b, theta_index_in_b])
      
      
      z_theta <- theta_val / (se_theta)
      
      
      p_theta <- 1 - pnorm(z_theta)
      
      
      
      cat(sprintf(" theta : %f (SE: %f) p = %f\n", theta_val, se_theta, p_theta))
      cat("\n")
    }
    
    
    # --- Weibull baseline hazard parameters ---
    cat("Scales and shapes of the Weibull baseline hazard\n")
    cat("------------\n")
    cat(sprintf(" %-8s %-12s %-12s %-12s %-12s\n", "Event", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
    
    current_b_idx_weib <- 1 
    for (k in 1:n_competing_events_from_formulas) {
      scale_val <- exp(fit1$b[current_b_idx_weib])
      se_scale <- sqrt(vcov_matrix[current_b_idx_weib, current_b_idx_weib])
      current_b_idx_weib <- current_b_idx_weib + 1
      
      shape_val <- exp(fit1$b[current_b_idx_weib])
      se_shape <- sqrt(vcov_matrix[current_b_idx_weib, current_b_idx_weib])
      current_b_idx_weib <- current_b_idx_weib + 1
      
      cat(sprintf(" 0->%-2d %12.5f %12.5f %12.5f %12.5f\n", k, scale_val, se_scale, shape_val, se_shape))
    }
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
    for (k in 1:n_competing_events_from_formulas) {
      cat(sprintf("     0 -> %-2d = %d\n", k, length(which(delta == k))))
    }
    cat("Lost to follow-up  = ", length(which(delta == 0)), "\n")
    cat("\n")
    
    if(fit1$istop==1){
      cat("Number of iterations: ", fit1$ni, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
      cat("All criteria were satisfied")
    }
    
    if(fit1$istop==2){
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
    
    
    if(length(partialH)>0){
      if(fit1$istop==3){
        cat("Number of iterations: ", fit1$ni, "\n")
        cat("Convergence criteria:\n")
        cat("------------\n")
        cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
        cat(sprintf("  %10.5f   %10.5f   %32.5f\n", fit1$ca, fit1$cb, fit1$rdm))
        cat("All criteria were satisfied, parameters in partialH were dropped from Hessian to define the relative distance to optimum")
      }
      
      
      if(fit1$istop==2){
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
  
  cat("\n")
  
  v <- fit1$v
  n <- nrow(data)
  p <- length(fit1$b)
  vcov_matrix <- matrix(0, p, p)
  vcov_matrix[upper.tri(vcov_matrix, diag = TRUE)] <- v
  
  vcov_matrix <- vcov_matrix + t(vcov_matrix) - diag(diag(vcov_matrix))
  
  
  jacobian_diag_elements <- c()
  current_b_idx <- 1 
  
  for (k in 1:n_competing_events_from_formulas) {
    jacobian_diag_elements <- c(jacobian_diag_elements, exp(fit1$b[current_b_idx]))
    current_b_idx <- current_b_idx + 1
    jacobian_diag_elements <- c(jacobian_diag_elements, exp(fit1$b[current_b_idx]))
    current_b_idx <- current_b_idx + 1
  }
  
  if (Frailty == 1) {
    jacobian_diag_elements <- c(jacobian_diag_elements, 2 * fit1$b[current_b_idx])
    current_b_idx <- current_b_idx + 1
  }
  
  total_beta_count <- sum(sapply(Xmat_list, ncol))
  if (total_beta_count > 0) {
    jacobian_diag_elements <- c(jacobian_diag_elements, rep(1, total_beta_count))
  }
  
  
  ## DELTA METHOD
  jacobian_matrix <- diag(jacobian_diag_elements)
  vcov_matrix <- jacobian_matrix %*% vcov_matrix %*% jacobian_matrix
  
  
  
  
  scale.weib <- numeric(n_competing_events_from_formulas)
  shape.weib <- numeric(n_competing_events_from_formulas)
  for (k in 1:n_competing_events_from_formulas) {
    scale.weib[k] <- exp(fit1$b[2*k -1])
    shape.weib[k] <- exp(fit1$b[2*k])
  }
  
  
  vcov <- vcov_matrix
  call <- original_call
  crit <- fit1$istop
  grad <- fit1$grad
  
  if(Frailty==1){
    groups <- unique(group)
  }
  
  
  
  n.events <- numeric(0)
  for (k in 1:n_competing_events_from_formulas) {
    n.events[paste0("0->", k)] <- length(which(delta == k))
  }
  n.events["Censored"] <- length(which(delta == 0))
  
  loglik <- fit1$fn.value
  
  
  
  
  
  beta_start_index_in_b <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1
  coef <- fit1$b[beta_start_index_in_b:p]
  names(coef) <- all_beta_names
  
  
  n.iter <- fit1$ni
  AIC <- (1/n) *(length(fit1$b) - fit1$fn.value)
  
  
  beta_p.value <- c()
  current_b_idx_for_betas_pvalue <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1
  current_vcov_idx_for_betas_pvalue <- length(weibull_start_vals) + (if (Frailty == 1) 1 else 0) + 1
  
  for (k in 1:n_competing_events_from_formulas) {
    current_Xmat <- Xmat_list[[k]]
    num_betas_k <- ncol(current_Xmat)
    
    if (num_betas_k > 0) {
      betas_k <- fit1$b[current_b_idx_for_betas_pvalue:(current_b_idx_for_betas_pvalue + num_betas_k - 1)]
      se_betas_k <- sqrt(diag(vcov_matrix)[current_vcov_idx_for_betas_pvalue:(current_vcov_idx_for_betas_pvalue + num_betas_k - 1)])
      z_values_k <- betas_k / se_betas_k
      p_values_k <- 2 * pnorm(abs(z_values_k), lower.tail = FALSE)
      beta_p.value <- c(beta_p.value, p_values_k)
    }
    current_b_idx_for_betas_pvalue <- current_b_idx_for_betas_pvalue + num_betas_k
    current_vcov_idx_for_betas_pvalue <- current_vcov_idx_for_betas_pvalue + num_betas_k
  }
  names(beta_p.value) <- all_beta_names 
  
  
  
  theta <- NULL
  VarTheta <- NULL
  theta_p.value <- NULL
  
  if (Frailty == 1) {
    theta_index_in_b <- length(weibull_start_vals) + 1
    theta <- (fit1$b[theta_index_in_b])^2
    VarTheta <- vcov_matrix[theta_index_in_b, theta_index_in_b] 
    
    se_theta_for_pvalue <- sqrt(VarTheta)
    z_theta_for_pvalue <- theta / (se_theta_for_pvalue)
    theta_p.value <- 1 - pnorm(z_theta_for_pvalue)
  }
  
  
  b <- c()
  for (k in 1:n_competing_events_from_formulas) {
    b <- c(b, scale.weib[k], shape.weib[k])
  }
  
  
  
  if(Frailty==1){
    b <- c(b,theta)
  }
  
  b <- c(b,coef)
  
  npar <- length(b)
  
  nvar=sum(sapply(Xmat_list, ncol))
  
  
  ##BASELINE HAZARDS AND CONFIDENCE BANDS IF POSSIBLE
  names(x0) <- paste0("x0.event", 1:n_competing_events_from_formulas)
  
  
  lam0 <- vector("list", n_competing_events_from_formulas)
  error_fit_haz <- vector("list", n_competing_events_from_formulas)
  
  for (k in 1:n_competing_events_from_formulas) {
    param_indices_k <- c(2*k - 1, 2*k) 
    
    if (any(param_indices_k %in% partialH)) {
      lam0[[k]] <- base_hazard(x0[[k]], scale.weib[k], shape.weib[k])
    } else {
      
      tryCatch({
        lam0[[k]] <- hazard_and_confidence_bands_cmprsk(x0[[k]], scale.weib[k], shape.weib[k], vcov_matrix[param_indices_k, param_indices_k])
      }, error = function(e) {
        error_fit_haz[[k]] <- e$message
      })
      
      if (!is.null(error_fit_haz[[k]])) {
        lam0[[k]] <- base_hazard(x0[[k]], scale.weib[k], shape.weib[k])
      }
    }
  }
  
  
  names(lam0) <- paste0("lam0.event", 1:n_competing_events_from_formulas)
  
  
  
  
  ####
  
  
  # --- Baseline Survivals and Confidence Bands ---
  surv0 <- vector("list", n_competing_events_from_formulas)
  
  medians <- c()
  lower_medians <- c()
  upper_medians <- c()
  
  error_fit_surv <- vector("list", n_competing_events_from_formulas)
  
  for (k in 1:n_competing_events_from_formulas) {
    param_indices_k <- c(2*k - 1, 2*k) 
    
    current_median <- NA
    current_lower <- NA
    current_upper <- NA
    
    if (any(param_indices_k %in% partialH)) {
      surv0[[k]] <- su(x0[[k]], scale.weib[k], shape.weib[k])
      current_median <- minmin(surv0[[k]][,2], x0[[k]])
    } else {
      tryCatch({
        surv0[[k]] <- survival_and_confidence_bands_cmprsk(x0[[k]], scale.weib[k], shape.weib[k], vcov_matrix[param_indices_k, param_indices_k])
      }, error = function(e) {
        error_fit_surv[[k]] <<- e$message
      })
      
      if (is.null(error_fit_surv[[k]])) {
        current_median <- minmin(surv0[[k]][,2], x0[[k]])
        current_lower <- minmin(surv0[[k]][,3], x0[[k]])
        current_upper <- minmin(surv0[[k]][,4], x0[[k]])
      } else {
        surv0[[k]] <- su(x0[[k]], scale.weib[k], shape.weib[k])
        current_median <- minmin(surv0[[k]][,2], x0[[k]])
        
      }
    }
    
    medians <- c(medians, current_median)
    lower_medians <- c(lower_medians, current_lower)
    upper_medians <- c(upper_medians, current_upper)
  }
  names(surv0) <- paste0("surv0.event", 1:n_competing_events_from_formulas)
  
  
  
  median_survivals_output <- vector("list", n_competing_events_from_formulas)
  for(k in 1:n_competing_events_from_formulas) {
    if (!is.na(lower_medians[k])) { 
      median_survivals_output[[k]] <- c(
        lower = lower_medians[k],
        median = medians[k],
        upper = upper_medians[k]
      )
    } else { 
      median_survivals_output[[k]] <- c(median = medians[k])
    }
    names(median_survivals_output)[k] <- paste0("event", k)
  }
  medians <- median_survivals_output
  
  
  
  ###
  
  ############### FRAILTY PREDICTION
  
  frailty.pred <- NULL
  frailty.var <- NULL
  frailty.sd <- NULL
  
  if (Frailty == 1) {
    eta <- vector("list", n_competing_events_from_formulas)
    q.y <- vector("list", n_competing_events_from_formulas)
    
    for (k in 1:n_competing_events_from_formulas) {
      current_Xmat <- Xmat_list[[k]]
      current_betas_k <- coef[grepl(paste0(".event", k), names(coef))] 
      
      eta[[k]] <- if (ncol(current_Xmat) == 0) rep(0, nrow(current_Xmat)) else as.vector(current_Xmat %*% current_betas_k)
      
      q.y[[k]] <- exp(eta[[k]]) * (y / scale.weib[k])^shape.weib[k]
      
      
      
      
      if(length(surv_trunc1)==3){ 
        q.y[[k]] <- q.y[[k]] - exp(eta[[k]]) * (l / scale.weib[k])^shape.weib[k]
      }
    }
    
    frailty.pred <- c()
    frailty.var <- c()
    frailty.sd <- c()
    
    unique_groups <- unique(group_data)
    for (j in unique_groups) {
      cluster_indices <- which(group_data == j)
      
      sum_delta_cluster <- 0
      for (k in 1:n_competing_events_from_formulas) {
        sum_delta_cluster <- sum_delta_cluster + sum(delta_indicators[[k]][cluster_indices])
      }
      
      sum_q.y_cluster <- 0
      for (k in 1:n_competing_events_from_formulas) {
        sum_q.y_cluster <- sum_q.y_cluster + sum(q.y[[k]][cluster_indices])
      }
      
      num.pred <- sum_delta_cluster + (1 / theta)
      denom.pred <- sum_q.y_cluster + (1 / theta)
      
      frailty.pred <- c(frailty.pred, num.pred / denom.pred)
      frailty.var <- c(frailty.var, num.pred / (denom.pred)^2)
      frailty.sd <- c(frailty.sd, sqrt(frailty.var[length(frailty.var)])) 
    }
    names(frailty.pred) <- unique_groups
    names(frailty.var) <- unique_groups
    names(frailty.sd) <- unique_groups
  }
  
  
  
  linear.pred <- vector("list", n_competing_events_from_formulas)
  
  for (k in 1:n_competing_events_from_formulas) {
    current_Xmat <- Xmat_list[[k]]
    current_betas_k <- coef[grepl(paste0(".event", k), names(coef))] 
    
    eta_k <- if (ncol(current_Xmat) == 0) rep(0, nrow(current_Xmat)) else as.vector(current_Xmat %*% current_betas_k)
    
    if (Frailty == 1) {
      individual_frailty_pred <- frailty.pred[match(group_data, names(frailty.pred))]
      linear.pred[[k]] <- eta_k + log(individual_frailty_pred)
    } else {
      linear.pred[[k]] <- eta_k
    }
    names(linear.pred[[k]]) <- rownames(data) 
  }
  names(linear.pred) <- paste0("linear.pred.event", 1:n_competing_events_from_formulas)
  
  
  
  
  
  
  Frailty <- ifelse(Frailty==1,TRUE,FALSE)
  
  
  
  
  
  
  
  
  
  
  ca= fit1$ca
  cb=fit1$cb
  rdm=fit1$rdm
  
  
  
  
  
  
  
  
  
  
  
  if(Frailty==TRUE){
    
    if(nvar>0){
      
      if(all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x0=x0,
                      lam0=lam0,
                      surv0=surv0,
                      medians=medians,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred=linear.pred
                      ,partialH=partialH,data=data_init
                      ,formulas=formulas
                      ,Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(!all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x0=x0,
                      lam0=lam0,
                      surv0=surv0, 
                      medians=medians, 
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred=linear.pred, 
                      partialH=partialH,
                      names.factor=names.factor,
                      global_chisq=global_chisq,
                      
                      p.global_chisq=p.global_chisq,
                      
                      dof_chisq=dof_chisq,
                      global_chisq.test=global_chisq.test,
                      data=data_init,
                      formulas=formulas,
                      Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    if(nvar==0){
      
      if(all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x0=x0,
                      lam0=lam0,
                      surv0=surv0,
                      medians=medians,
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred=linear.pred
                      ,partialH=partialH,data=data_init
                      ,formulas=formulas
                      ,Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(!all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,groups=groups,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,theta=theta,
                      VarTheta=VarTheta,Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,theta_p.value=theta_p.value,x0=x0,
                      lam0=lam0,
                      surv0=surv0, 
                      medians=medians, 
                      frailty.pred=frailty.pred,frailty.var=frailty.var,frailty.sd=frailty.sd,
                      linear.pred=linear.pred, 
                      partialH=partialH,
                      names.factor=names.factor,
                      global_chisq=global_chisq,
                      
                      p.global_chisq=p.global_chisq,
                      
                      dof_chisq=dof_chisq,
                      global_chisq.test=global_chisq.test,
                      data=data_init,
                      formulas=formulas,
                      Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    
    
  }
  
  if(Frailty==FALSE){
    
    if(nvar>0){
      
      if(all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,
                      Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x0=x0,
                      lam0=lam0,
                      surv0=surv0,
                      medians=medians,
                      
                      linear.pred=linear.pred
                      ,partialH=partialH,data=data_init
                      ,formulas=formulas
                      ,Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(!all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,coef=coef,
                      n.iter=n.iter,AIC=AIC,beta_p.value=beta_p.value,
                      Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x0=x0,
                      lam0=lam0,
                      surv0=surv0, 
                      medians=medians, 
                      
                      linear.pred=linear.pred, 
                      partialH=partialH,
                      names.factor=names.factor,
                      global_chisq=global_chisq,
                      
                      p.global_chisq=p.global_chisq,
                      
                      dof_chisq=dof_chisq,
                      global_chisq.test=global_chisq.test,
                      data=data_init,
                      formulas=formulas,
                      Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    if(nvar==0){
      
      if(all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,
                      Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x0=x0,
                      lam0=lam0,
                      surv0=surv0,
                      medians=medians,
                      
                      linear.pred=linear.pred
                      ,partialH=partialH,data=data_init
                      ,formulas=formulas
                      ,Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
      
      
      if(!all(sapply(factor_covariates_per_event, is.null))){
        value <- list(b=b,scale.weib=scale.weib,shape.weib=shape.weib,vcov=vcov,call=call,
                      n=n,n.events=n.events,loglik=loglik,
                      n.iter=n.iter,AIC=AIC,
                      Frailty=Frailty,crit=crit,grad=grad,npar=npar,
                      nvar=nvar,x0=x0,
                      lam0=lam0,
                      surv0=surv0, 
                      medians=medians, 
                      
                      linear.pred=linear.pred, 
                      partialH=partialH,
                      names.factor=names.factor,
                      global_chisq=global_chisq,
                      
                      p.global_chisq=p.global_chisq,
                      
                      dof_chisq=dof_chisq,
                      global_chisq.test=global_chisq.test,
                      data=data_init,
                      formulas=formulas,
                      Xmat=Xmat_list,
                      delta=delta,trunc=trunc,
                      ca=ca,cb=cb,rdm=rdm,LIMparam=LIMparam,LIMlogl=LIMlogl,LIMderiv=LIMderiv)
      }
    }
    
    
    
    
  }
  
  
  class(value)="frailtyCmprsk"
  
  return(invisible(value))}








