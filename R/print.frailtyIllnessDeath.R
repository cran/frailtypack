#' Print a Short Summary of parameter estimates of a Weibull Illness-Death model with (or without) shared frailty between transitions.
#' 
#' Prints a short summary of parameter estimates of a 'frailtyIllnessDeath' object
#' 
#' 
#' @usage \method{print}{frailtyIllnessDeath}(x,
#' ...)
#' @param x the result of a call to the frailtyIllnessDeath function.
#' @param \dots other unused arguments.
#' @return
#' 
#' Print the parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{frailtyIllnessDeath}}
#' @keywords methods
##' @export
"print.frailtyIllnessDeath" <- function (x, ...) 
{
  
  
  
  

  if (!inherits(x, "frailtyIllnessDeath")){
    stop("Object must be of class 'frailtyIllnessDeath'")
  }

 
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
  
  
if(x$model == "Semi-Markov")
{
  
  
  
  
  if(any(grepl("cluster\\(.*\\)", deparse(x$formula)))){
    
    cat("Illness-death model with shared frailty between transitions\n")
    cat("Using Weibull baseline hazard functions \n")
    
    
    if( x$trunc==TRUE){
      cat("Left truncation structure is used\n")
    }
    
    cat("Semi-Markov model is used for the transition 1->2\n")
    cat("\n")
  }
  
  if(!any(grepl("cluster\\(.*\\)", deparse(x$formula)))){
    
    cat("Illness-death model\n")
    cat("Using Weibull baseline hazard functions \n")
    
    
    if( x$trunc==TRUE){
      cat("Left truncation structure is used\n")
    }
    
    cat("Semi-Markov model is used for the transition 1->2\n")
    cat("\n")
  }
  
}













if(x$model == "Markov")
{
  
  
  
  
  if(any(grepl("cluster\\(.*\\)", deparse(x$formula)))){
    
    cat("Illness-death model with shared frailty between transitions\n")
    cat("Using Weibull baseline hazard functions \n")
    
    
    if(x$trunc==TRUE){
      cat("Left truncation structure is used\n")
    }
    
    cat("Markov model is used for the transition 1->2\n")
    cat("\n")
  }
  
  
  if(!any(grepl("cluster\\(.*\\)", deparse(x$formula)))){
    
    cat("Illness-death model\n")
    cat("Using Weibull baseline hazard functions \n")
    
    
    if(x$trunc==TRUE){
      cat("Left truncation structure is used\n")
    }
    
    cat("Markov model is used for the transition 1->2\n")
    cat("\n")
  }
}





  if(x$crit==1){
   
    

   vcov_matrix <- x$vcov
    terms_object <- terms(x$formula)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_01 <- ifelse(x$Frailty == TRUE, 8, 7)
    
    if(length(covariates_without_cluster01)>0){
      covariates_transitions_01 <- x$b[debut_reg_01:(debut_reg_01 + ncol(x$Xmat1)-1)]
      
      cat("Transition 0 -> 1:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat1)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat1))){
        
        cov_name <- colnames(x$Xmat1)[i]
        
        
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
      
      if (length(x$global_chisq.01) > 0) {
        # Filter factors to only include those with dof > 1, as per standard global test definition
        factors_with_mult_dof <- names(x$dof_chisq.01)[x$dof_chisq.01 > 1]
        
        if (length(factors_with_mult_dof) > 0) { # Only print header if there are such factors
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.01[factor_name],
                        x$dof_chisq.01[factor_name],
                        x$p.global_chisq.01[factor_name]
            ))
          }
          cat("\n") 
        }
      }
      
    }
    
    
    
    
    
    
    
    terms_object <- terms(x$formula.terminalEvent)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_02 <- debut_reg_01 + ncol(x$Xmat1)
    if(length(covariates_without_cluster02)>0){
      
      covariates_transitions_02 <- x$b[debut_reg_02:(debut_reg_02 + ncol(x$Xmat2)-1)]
      
      cat("Transition 0 -> 2:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat2)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat2))){
        
        cov_name <- colnames(x$Xmat2)[i]
        
        
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
      
      if (length(x$global_chisq.02) > 0) {
        
        factors_with_mult_dof <- names(x$dof_chisq.02)[x$dof_chisq.02 > 1]
        
        if (length(factors_with_mult_dof) > 0) { 
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.02[factor_name],
                        x$dof_chisq.02[factor_name],
                        x$p.global_chisq.02[factor_name]
            ))
          }
          cat("\n") 
        }
      }
    }
    
    
    
    
    
    
    terms_object <- terms(x$formula.terminalEvent)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_12 <- debut_reg_02 + ncol(x$Xmat2)
    if(length(covariates_without_cluster12)>0){
      
      covariates_transitions_12 <- x$b[debut_reg_12:(debut_reg_12 + ncol(x$Xmat3)-1)]
      
      cat("Transition 1 -> 2:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat3)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat3))){
        
        cov_name <- colnames(x$Xmat3)[i]
        
        
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
      
      if (length(x$global_chisq.12) > 0) {
        
        factors_with_mult_dof <- names(x$dof_chisq.12)[x$dof_chisq.12 > 1]
        
        if (length(factors_with_mult_dof) > 0) { 
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.12[factor_name],
                        x$dof_chisq.12[factor_name],
                        x$p.global_chisq.12[factor_name]
            ))
          }
          cat("\n") 
        }
      }
    }
    
    
    
    if(x$Frailty == TRUE){
      cat("Frailty Parameters:\n")
      cat("------------\n")
      cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", x$b[7], sqrt(vcov_matrix[7,7]), 
                  1- pnorm((x$b[7]) / sqrt(vcov_matrix[7,7]))))
      cat("\n")
    }
    
    cat("Scales and shapes of the Weibull baseline hazard\n")
    cat("------------\n")
    cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
    
    
    cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[1], sqrt(vcov_matrix[1,1]),
                x$b[2], sqrt(vcov_matrix[2,2])))
    cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[3], sqrt(vcov_matrix[3,3]),
                x$b[4], sqrt(vcov_matrix[4,4])))
    cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[5], sqrt(vcov_matrix[5,5]),
                x$b[6], sqrt(vcov_matrix[6,6])))
    
    cat("\n")
    
    cat("The expression of the Weibull hazard function is: \n")
    cat("           'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)' \n" )
    cat("The expression of the Weibull survival function is: \n")
    cat("           'S(t) = exp[- (t/scale)^shape]' \n" )
    
    cat("\n")
    
    cat("Marginal log-likelihood = ", x$loglik, "\n")
    cat("AIC = Akaike information Criterion = ", (1/nrow(x$data)) * (length(x$b) - x$loglik), "\n")
    cat("           'AIC = (1/n)[np - l(.)]' \n")
    
    cat("\n")
    
    cat("Number of subjects = ", nrow(x$data), "\n")
    cat("             0 -> 1 = ", length(which(x$delta1 == 1)), "\n")
    cat("             0 -> 2 = ", length(which(x$delta1 == 0 & x$delta2 == 1)), "\n")
    cat("             1 -> 2 = ", length(which(x$delta1 == 1 & x$delta2 == 1)), "\n")
    cat("Lost to follow-up  = ", length(which(x$delta1 == 0 & x$delta2 == 0)), "\n")
    cat("\n")
    
    cat("Number of iterations: ", x$n.iter, "\n")
    cat("Convergence criteria:\n")
    cat("------------\n")
    cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
    cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
    cat("All criteria were satisfied")
    
  }


  if(x$crit==2){
   
   vcov_matrix <- x$vcov
    
    
    
    terms_object <- terms(x$formula)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_01 <- ifelse(x$Frailty == TRUE, 8, 7)
    if(length(covariates_without_cluster01)>0){
      
      covariates_transitions_01 <- x$b[debut_reg_01:(debut_reg_01 + ncol(x$Xmat1)-1)]
      
      cat("Transition 0 -> 1:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat1)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat1))){
        
        cov_name <- colnames(x$Xmat1)[i]
        
        
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
      
      if (length(x$global_chisq.01) > 0) {
        
        factors_with_mult_dof <- names(x$dof_chisq.01)[x$dof_chisq.01 > 1]
        
        if (length(factors_with_mult_dof) > 0) { 
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.01[factor_name],
                        x$dof_chisq.01[factor_name],
                        x$p.global_chisq.01[factor_name]
            ))
          }
          cat("\n") 
        }
      }
      
    }
    
    
    
    
    
    
    terms_object <- terms(x$formula.terminalEvent)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_02 <- debut_reg_01 + ncol(x$Xmat1)
    if(length(covariates_without_cluster02)>0){
      
      covariates_transitions_02 <- x$b[debut_reg_02:(debut_reg_02 + ncol(x$Xmat2)-1)]
      
      cat("Transition 0 -> 2:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat2)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat2))){
        
        cov_name <- colnames(x$Xmat2)[i]
        
        
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
      
      
      if (length(x$global_chisq.02) > 0) {
        
        factors_with_mult_dof <- names(x$dof_chisq.02)[x$dof_chisq.02 > 1]
        
        if (length(factors_with_mult_dof) > 0) { 
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.02[factor_name],
                        x$dof_chisq.02[factor_name],
                        x$p.global_chisq.02[factor_name]
            ))
          }
          cat("\n") 
        }
      }
      
    }
    
    
    
    
    
    
    terms_object <- terms(x$formula.terminalEvent)
    covariates <- attr(terms_object, "term.labels")
    covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
    debut_reg_12 <- debut_reg_02 + ncol(x$Xmat2)
    if(length(covariates_without_cluster12)>0){
      covariates_transitions_12 <- x$b[debut_reg_12:(debut_reg_12 + ncol(x$Xmat3)-1)]
      
      cat("Transition 1 -> 2:\n")
      cat("------------\n")
      
      max_cov_name_length <- max(nchar(colnames(x$Xmat3)))
      
      
      cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                  max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
      
      
      for(i in 1:length(colnames(x$Xmat3))){
        
        cov_name <- colnames(x$Xmat3)[i]
        
        
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
      
      if (length(x$global_chisq.12) > 0) {
        
        factors_with_mult_dof <- names(x$dof_chisq.12)[x$dof_chisq.12 > 1]
        
        if (length(factors_with_mult_dof) > 0) { 
          cat(sprintf(" %-*s %12s %12s %12s \n",
                      max_cov_name_length, "", "chisq", "df", "global p"))
          
          for (factor_name in factors_with_mult_dof) {
            cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                        max_cov_name_length, factor_name,
                        x$global_chisq.12[factor_name],
                        x$dof_chisq.12[factor_name],
                        x$p.global_chisq.12[factor_name]
            ))
          }
          cat("\n") 
        }
      }
      
    }
    
    
    
    # Frailty Parameters
    if(x$Frailty == TRUE){
      cat("Frailty Parameters:\n")
      cat("------------\n")
      cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", x$b[7], sqrt(vcov_matrix[7,7]), 
                  1- pnorm((x$b[7]) / sqrt(vcov_matrix[7,7]))))
      cat("\n")
    }
    
    # Weibull baseline hazard parameters
    cat("Scales and shapes of the Weibull baseline hazard\n")
    cat("------------\n")
    cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
    
    
    cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[1], sqrt(vcov_matrix[1,1]),
                x$b[2], sqrt(vcov_matrix[2,2])))
    cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[3], sqrt(vcov_matrix[3,3]),
                x$b[4], sqrt(vcov_matrix[4,4])))
    cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[5], sqrt(vcov_matrix[5,5]),
                x$b[6], sqrt(vcov_matrix[6,6])))
    
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
    cat("             0 -> 1 = ", length(which(x$delta1 == 1)), "\n")
    cat("             0 -> 2 = ", length(which(x$delta1 == 0 & x$delta2 == 1)), "\n")
    cat("             1 -> 2 = ", length(which(x$delta1 == 1 & x$delta2 == 1)), "\n")
    cat("Lost to follow-up  = ", length(which(x$delta1 == 0 & x$delta2 == 0)), "\n")
    cat("\n")
    
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
      vcov_matrix <- x$vcov
      
      terms_object <- terms(x$formula)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_01 <- ifelse(x$Frailty == TRUE, 8, 7)
      if(length(covariates_without_cluster01)>0){
        
        covariates_transitions_01 <- x$b[debut_reg_01:(debut_reg_01 + ncol(x$Xmat1)-1)]
        
        cat("Transition 0 -> 1:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat1)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat1))){
          
          cov_name <- colnames(x$Xmat1)[i]
          
          
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
        
        if (length(x$global_chisq.01) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.01)[x$dof_chisq.01 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.01[factor_name],
                          x$dof_chisq.01[factor_name],
                          x$p.global_chisq.01[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      
      
      terms_object <- terms(x$formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_02 <- debut_reg_01 + ncol(x$Xmat1)
      if(length(covariates_without_cluster02)>0){
        
        covariates_transitions_02 <- x$b[debut_reg_02:(debut_reg_02 + ncol(x$Xmat2)-1)]
        
        cat("Transition 0 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat2)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat2))){
          
          cov_name <- colnames(x$Xmat2)[i]
          
          
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
        
        
        if (length(x$global_chisq.02) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.02)[x$dof_chisq.02 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.02[factor_name],
                          x$dof_chisq.02[factor_name],
                          x$p.global_chisq.02[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
        
      }
      
      
      
      
      
      terms_object <- terms(x$formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_12 <- debut_reg_02 + ncol(x$Xmat2)
      if(length(covariates_without_cluster12)>0){
        
        covariates_transitions_12 <- x$b[debut_reg_12:(debut_reg_12 + ncol(x$Xmat3)-1)]
        
        cat("Transition 1 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat3)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat3))){
          
          cov_name <- colnames(x$Xmat3)[i]
          
          
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
        
        
        if (length(x$global_chisq.12) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.12)[x$dof_chisq.12 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.12[factor_name],
                          x$dof_chisq.12[factor_name],
                          x$p.global_chisq.12[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      # Frailty Parameters
      if(x$Frailty == TRUE){
        cat("Frailty Parameters:\n")
        cat("------------\n")
        cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", x$b[7], sqrt(vcov_matrix[7,7]), 
                    1- pnorm((x$b[7]) / sqrt(vcov_matrix[7,7]))))
        cat("\n")
      }
      
      # Weibull baseline hazard parameters
      cat("Scales and shapes of the Weibull baseline hazard\n")
      cat("------------\n")
      cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
      
      
      cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[1], sqrt(vcov_matrix[1,1]),
                  x$b[2], sqrt(vcov_matrix[2,2])))
      cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[3], sqrt(vcov_matrix[3,3]),
                  x$b[4], sqrt(vcov_matrix[4,4])))
      cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[5], sqrt(vcov_matrix[5,5]),
                  x$b[6], sqrt(vcov_matrix[6,6])))
      
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
      cat("             0 -> 1 = ", length(which(x$delta1 == 1)), "\n")
      cat("             0 -> 2 = ", length(which(x$delta1 == 0 & x$delta2 == 1)), "\n")
      cat("             1 -> 2 = ", length(which(x$delta1 == 1 & x$delta2 == 1)), "\n")
      cat("Lost to follow-up  = ", length(which(x$delta1 == 0 & x$delta2 == 0)), "\n")
      cat("\n")
      
      cat("Number of iterations: ", x$n.iter, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
      cat("All criteria were satisfied, parameters in partialH were dropped from Hessian to define the relative distance to optimum")
      
    }
    
  }



  if(length(x$partialH)>0){
    if(x$crit==2){
     vcov_matrix <- x$vcov
      
      terms_object <- terms(x$formula)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster01 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_01 <- ifelse(x$Frailty == TRUE, 8, 7)
      
      if(length(covariates_without_cluster01)>0){
        
        covariates_transitions_01 <- x$b[debut_reg_01:(debut_reg_01 + ncol(x$Xmat1)-1)]
        
        cat("Transition 0 -> 1:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat1)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat1))){
          
          cov_name <- colnames(x$Xmat1)[i]
          
          
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
        
        if (length(x$global_chisq.01) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.01)[x$dof_chisq.01 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.01[factor_name],
                          x$dof_chisq.01[factor_name],
                          x$p.global_chisq.01[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      
      terms_object <- terms(x$formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster02 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_02 <- debut_reg_01 + ncol(x$Xmat1)
      if(length(covariates_without_cluster02)>0){
        
        covariates_transitions_02 <- x$b[debut_reg_02:(debut_reg_02 +ncol(x$Xmat2)-1)]
        
        cat("Transition 0 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat2)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat2))){
          
          cov_name <- colnames(x$Xmat2)[i]
          
          
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
        
        if (length(x$global_chisq.02) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.02)[x$dof_chisq.02 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.02[factor_name],
                          x$dof_chisq.02[factor_name],
                          x$p.global_chisq.02[factor_name]
              ))
            }
            cat("\n") 
          }
        }
        
      }
      
      
      
      
      terms_object <- terms(x$formula.terminalEvent)
      covariates <- attr(terms_object, "term.labels")
      covariates_without_cluster12 <- covariates[!grepl("cluster\\(.*\\)", covariates)]
      debut_reg_12 <- debut_reg_02 + ncol(x$Xmat2)
      if(length(covariates_without_cluster12)>0){
        
        covariates_transitions_12 <- x$b[debut_reg_12:(debut_reg_12 + ncol(x$Xmat3)-1)]
        
        cat("Transition 1 -> 2:\n")
        cat("------------\n")
        
        max_cov_name_length <- max(nchar(colnames(x$Xmat3)))
        
        
        cat(sprintf(" %-*s %12s %12s %12s %12s %12s\n", 
                    max_cov_name_length, "", "coef", "exp(coef)", "SE(coef)", "z", "p"))
        
        
        for(i in 1:length(colnames(x$Xmat3))){
          
          cov_name <- colnames(x$Xmat3)[i]
          
          
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
        
        if (length(x$global_chisq.12) > 0) {
          
          factors_with_mult_dof <- names(x$dof_chisq.12)[x$dof_chisq.12 > 1]
          
          if (length(factors_with_mult_dof) > 0) { 
            cat(sprintf(" %-*s %12s %12s %12s \n",
                        max_cov_name_length, "", "chisq", "df", "global p"))
            
            for (factor_name in factors_with_mult_dof) {
              cat(sprintf(" %-*s %12.5f %12.5f %12.5f\n",
                          max_cov_name_length, factor_name,
                          x$global_chisq.12[factor_name],
                          x$dof_chisq.12[factor_name],
                          x$p.global_chisq.12[factor_name]
              ))
            }
            cat("\n") 
          }
        }
      }
      
      
      # Frailty Parameters
      if(x$Frailty == TRUE){
        cat("Frailty Parameters:\n")
        cat("------------\n")
        cat(sprintf(" theta : %f (SE (H): %f) p = %f\n", x$b[7], sqrt(vcov_matrix[7,7]), 
                    1- pnorm((x$b[7]) / sqrt(vcov_matrix[7,7]))))
        cat("\n")
      }
      
      # Weibull baseline hazard parameters
      cat("Scales and shapes of the Weibull baseline hazard\n")
      cat("------------\n")
      cat(sprintf("            %-8.5s %-16s %-8.5s %-10s\n", "Scale", "SE(Scale)", "Shape", "SE(Shape)"))
      
      
      cat(sprintf(" 0->1   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[1], sqrt(vcov_matrix[1,1]),
                  x$b[2], sqrt(vcov_matrix[2,2])))
      cat(sprintf(" 0->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[3], sqrt(vcov_matrix[3,3]),
                  x$b[4], sqrt(vcov_matrix[4,4])))
      cat(sprintf(" 1->2   %10.5f   %10.5f   %10.5f   %10.5f\n", x$b[5], sqrt(vcov_matrix[5,5]),
                  x$b[6], sqrt(vcov_matrix[6,6])))
      
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
      cat("             0 -> 1 = ", length(which(x$delta1 == 1)), "\n")
      cat("             0 -> 2 = ", length(which(x$delta1 == 0 & x$delta2 == 1)), "\n")
      cat("             1 -> 2 = ", length(which(x$delta1 == 1 & x$delta2 == 1)), "\n")
      cat("Lost to follow-up  = ", length(which(x$delta1 == 0 & x$delta2 == 0)), "\n")
      cat("\n")
      
      cat("Number of iterations: ", x$n.iter, "\n")
      cat("Convergence criteria:\n")
      cat("------------\n")
      cat(sprintf("  %-14s %-14s %-7s\n" ,"Parameters","Function","Relative Distance to optimum"))
      cat(sprintf("  %10.5f   %10.5f   %32.5f\n", x$ca, x$cb, x$rdm))
      cat("Maximum number of iterations reached, parameters in partialH were dropped from Hessian to define the relative distance to optimum\n")
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


}
