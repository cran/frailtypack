#' Print a Short Summary of parameter estimates of a Weibull competing risks model with (or without) shared frailty between transitions.
#' 
#' Prints a short summary of parameter estimates of a 'frailtyCmprsk' object
#' 
#' 
#' @usage \method{print}{frailtyCmprsk}(x,
#' ...)
#' @param x the result of a call to the frailtyCmprsk function.
#' @param \dots other unused arguments.
#' @return
#' 
#' Print the parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{frailtyCmprsk}}
#' @keywords methods
##' @export
"print.frailtyCmprsk" <- function (x, ...) 
{
  if (!inherits(x, "frailtyCmprsk")) 
    stop("Object must be of class 'frailtyCmprsk'")
  
  
  
  
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  
  if(x$Frailty==TRUE){
    
    cat("Competing risks model with shared frailty between transitions\n")
    cat("Using Weibull baseline hazard functions \n")
    if(x$trunc==TRUE ){
      cat("Left truncation structure is used\n")
    }
    
    
    cat("\n")
  }
  
  if(x$Frailty==FALSE){
    
    cat("Competing risks model\n")
    cat("Using Weibull baseline hazard functions \n")
    if(x$trunc==TRUE){
      cat("Left truncation structure is used\n")
    }
    
    
    cat("\n")
    
  }
  
  vcov_matrix <- x$vcov
  
  
   n_competing_events_from_formulas <- length(x$formulas)
  

  weibull_start_vals <- c(x$scale.weib,x$shape.weib)

  beta_start_index_in_b <- length(weibull_start_vals) + (if (x$Frailty == TRUE) 1 else 0) + 1
  
  current_beta_idx_in_b <- beta_start_index_in_b
  current_vcov_idx_for_factors <- beta_start_index_in_b
  
  
    
    
    
    current_b_idx_for_betas <- length(weibull_start_vals) + (if (x$Frailty == TRUE) 1 else 0) + 1 
    current_vcov_idx_for_betas <- length(weibull_start_vals) + (if (x$Frailty == TRUE) 1 else 0) + 1 
    
    for (k in 1:n_competing_events_from_formulas) {
      current_Xmat <- x$Xmat[[k]]
      num_betas_k <- ncol(current_Xmat)
      
      if (num_betas_k > 0) {
        cat(paste0("Transition 0 -> ", k, ":\n"))
        cat("------------\n")
        
        betas_k <- x$b[current_b_idx_for_betas:(current_b_idx_for_betas + num_betas_k - 1)]
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
        if (length(x$global_chisq[[k]]) > 0) {
          factors_with_mult_dof <- names(x$dof_chisq[[k]])[x$dof_chisq[[k]] > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq[[k]][factor_name],
                          x$dof_chisq[[k]][factor_name],
                          x$p.global_chisq[[k]][factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
        current_b_idx_for_betas <- current_b_idx_for_betas + num_betas_k
        current_vcov_idx_for_betas <- current_vcov_idx_for_betas + num_betas_k
      }
    }
    
    
    if (x$Frailty == TRUE) {
      cat("Frailty Parameters:\n")
      cat("------------\n")
      theta_index_in_b <- length(weibull_start_vals) + 1 
      
      theta_val <- x$b[theta_index_in_b]
      
      
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
      scale_val <- x$b[current_b_idx_weib]
      se_scale <- sqrt(vcov_matrix[current_b_idx_weib, current_b_idx_weib])
      current_b_idx_weib <- current_b_idx_weib + 1
      
      shape_val <- x$b[current_b_idx_weib]
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
    cat("Marginal log-likelihood = ", x$loglik, "\n")
    cat("AIC = Akaike information Criterion = ", (1/nrow(x$data)) * (length(x$b) - x$loglik), "\n")
    cat("           'AIC = (1/n)[np - l(.)]' \n")
    
    cat("\n")
    
    # Number of subjects and transitions
    cat("Number of subjects = ", nrow(x$data), "\n")
    for (k in 1:n_competing_events_from_formulas) {
      cat(sprintf("     0 -> %-2d = %d\n", k, length(which(x$delta == k))))
    }
    cat("Lost to follow-up  = ", length(which(x$delta == 0)), "\n")
    cat("\n")
    
    if(x$crit==1){
      cat("Number of iterations: ", x$n.iter, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
      cat("All criteria were satisfied")
    }
    
    if(x$crit==2){
      cat("Number of iterations: ", x$n.iter, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
      cat("Maximum number of iterations reached\n")
      if (x$ca > x$LIMparam) {
        cat("  Parameter criterion not reached. Threshold: ", x$LIMparam, "\n")
      }
      
      if (x$cb > x$LIMlogl) {
        cat("  Function criterion not reached. Threshold: ", x$LIMlogl, "\n")
      }
      
      if (x$rdm > x$LIMderiv) {
        cat("  Relative distance to optimum criterion not reached. Threshold: ", x$LIMderiv, "\n")
      }
    }
    
    
    if(length(x$partialH)>0){
      if(x$crit==3){
        cat("Number of iterations: ", x$n.iter, "\n")
        cat("Convergence criteria:\n")
        cat("------------\n")
        cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
        cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
        cat("All criteria were satisfied, parameters in x$partialH were dropped from Hessian to define the relative distance to optimum")
      }
      
      
      if(x$crit==2){
        cat("Number of iterations: ", x$n.iter, "\n")
        cat("Convergence criteria:\n")
        cat("------------\n")
        cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
        cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
        cat("Maximum number of iterations reached, parameters in x$partialH were dropped from Hessian to define the relative distance to optimum\n")
        if (x$ca > x$LIMparam) {
          cat("  Parameter criterion not reached. Threshold: ", x$LIMparam, "\n")
        }
        
        if (x$cb > x$LIMlogl) {
          cat("  Function criterion not reached. Threshold: ", x$LIMlogl, "\n")
        }
        
        if (x$rdm > x$LIMderiv) {
          cat("  Relative distance to optimum criterion not reached. Threshold: ", x$LIMderiv, "\n")
        }
      }
      
    }
  
  
  cat("\n")
  
  
}

