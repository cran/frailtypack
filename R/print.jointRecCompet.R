#' Print a Short Summary of parameter estimates of a joint competing risks midel
#' 
#' Prints a short summary of parameter estimates of a joint competing risks model
#' or more generally an object of class 'jointRecCompet'.
#' 
#' 
#' @usage \method{print}{jointRecCompet}(x, digits = max(options()$digits - 4, 6),
#' ...)
#' @param x the result of a call to the jointRecCompet function
#' @param digits number of digits to print
#' @param \dots other unused arguments
#' @return
#' 
#' Print, separately for each type of event (Recurrent, Terminal1 and
#' Terminal2), the parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{jointRecCompet}}
#' @keywords methods jointRecCompet
#' @export
"print.jointRecCompet" <-function(x, digits = max(options()$digits - 4, 6), ...){
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$summary.table$Estimate
  nvar <- length(coef)
  if (x$critCV[2] == 1){
    if (!is.null(coef)){
      
      seH <- sqrt(diag(x$varH.Estimate))
      cat("\n")
      cat(" Joint competing risk model for one recurrent event and two terminal events \n")
      if(x$controls[3] == 0){
        if(x$controls[4]==1) cat(" Using splines baseline hazard functions using equidistant intervals \n")
        if(x$controls[4]==0) cat(" Using splines baseline hazard functions using percentile-based intervals \n")
      }else{
        cat(" Using Weibull baseline hazard functions","\n")
      }
      if (!is.null(cl <- x$call)){
        cat("\n")
        
        cat("Call:\n")
        dput(cl)
        
        cat("\n")
      }
      
      if(x$noVarEvent[1] == 0){
        
        cat("\n")
        cat("Recurrent event:\n")
        cat("------------ \n")
        xx=grep("Recurrent", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Recurrent: ","",tmp$Parameter)
        tmp$H0=gsub("Recurrent: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          xx=grep("Weibull", tmp[,1])
          print(tmp[-xx,], , drop=FALSE,row.names=FALSE)
        }else{
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          print(tmp, , drop=FALSE,row.names=FALSE)
        }
                
      }
      
      
      if(x$noVarEvent[2] == 0){	
        cat("Terminal event 1:\n")
        cat("------------- \n")
        xx=grep("Terminal1", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Terminal1: ","",tmp$Parameter)
        tmp$H0=gsub("Terminal1: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          xx1=grep("Weibull", tmp[,1])
          xx2 = grep("Alpha", tmp[,1])
          print(tmp[-c(xx1,xx2),], , drop=FALSE,row.names=FALSE)
        }else{
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          xx = grep("Alpha", tmp[,1])
          print(tmp[-xx,], , drop=FALSE,row.names=FALSE)
        }

        cat("\n")
      }
      
      
      if(x$noVarEvent[4]== 0){	
        cat("Terminal event 2:\n")
        cat("------------- \n")
        xx=grep("Terminal2", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Terminal2: ","",tmp$Parameter)
        tmp$H0=gsub("Terminal2: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          xx1=grep("Weibull", tmp[,1])
          xx2 = grep("Alpha", tmp[,1])
          print(tmp[-c(xx1,xx2),], , drop=FALSE,row.names=FALSE)
        }else{
          colnames(tmp) = c("","coef","exp(coef)","SE coef","z","p")
          xx = grep("Alpha", tmp[,1])
          print(tmp[-xx,], , drop=FALSE,row.names=FALSE)
        }

        cat("\n")
      }
      
      
    }
    
    #sigma
    tmp = x$summary.table
    xsig=grep("Sigma",tmp$Parameter)
    sigma <- tmp[xsig,"Estimate"]
    sesigma <- tmp[xsig,"Estimate.SE"]
    
    
    #AD:
    if (x$noVarEvent[1] == 1){
      cat("\n")
      cat("    Recurrent event: No covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVarEvent[2] == 1){
      cat("\n")
      cat("    Terminal event 1: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    if (x$noVarEvent[4] == 1){
      cat("\n")
      cat("    Terminal event 2: No covariates \n")
      cat("    ----------- \n")
    }
    #AD:  
    cat("Frailty parameters: \n")
    cat("----------- \n")
    
    cat(" theta :", sigma, "(SE (H):", sesigma, ")", "p =", ifelse(signif(1 - pnorm(sigma/sesigma), digits - 1) == 0, "< 1e-16", signif(1 - pnorm(sigma/sesigma))), "\n")
    cat(" \n")
    cat("Power parameters (alphas)","\n")
    
    tmpt1 = as.numeric(tmp[grep("Terminal1: Alpha",tmp[,1]),c(4,5,8)])
    tmpt2 = as.numeric(tmp[grep("Terminal2: Alpha",tmp[,1]),c(4,5,8)])
    cat("    Terminal event 1: ", tmpt1[1],"(SE: ",tmpt1[2],")","p = ",tmpt1[3])
    cat('\n')
    cat("    Terminal event 2: ", tmpt2[1],"(SE: ",tmpt2[2],")","p = ",tmpt2[3])
    cat('\n') 
    
    
    
    
    
    if(x$controls[3] == 0){
      cat("\n")
      cat(paste("Penalized marginal log-likelihood =", round(x$logLikPenal,2)))
      cat("\n")
      cat("LCV = the approximate likelihood cross-validation criterion\n")
      cat("  in the semi parametric case     =",x$LCV,"\n")
    }else{
      cat("\n")
      
      cat("Scales and shapes of the Weibull baseline hazard","\n")
      
      xx=which(tmp$Parameter %in% c("Recurrent: Scale","Recurrent: Shape",
                                    "Terminal1: Scale", "Terminal1: Shape",
                                    "Terminal2: Scale","Terminal2: Shape"))
      tmphz=data.frame(matrix(tmp[xx,"Estimate"],3,2,byrow = T))
      
      rownames(tmphz)=c(" Recurrent", " Terminal 1"," Terminal 2")
      colnames(tmphz) =c("Scale","Shape")
      
      print(tmphz, , drop=FALSE)
      
      cat("\n")
      cat(paste("Marginal log-likelihood =", round(x$logLik,2)))
      cat("\n")
      cat(paste("Akaike information Criterion  =", round(x$AIC,2)))
      cat("\n")
      
    }
    
    #AD:
    cat("\n")
    cat("Number of subjects = ", x$n) 
    cat("\n")
    cat("Number of recurrent events = ", x$nevts[1])
    cat("\n")
    cat("Number of terminal events of type 1 = ", x$nevts[2])
    cat("\n")
    cat("Number of terminal events of type 2 = ", x$nevts[3])
    cat("\n")
    cat("Number of iterations: ", x$ni,"\n")
    cat("Number of nodes for the Gauss-Herminte quadrature:",length(x$ghNodes))
    if(x$controls[3] == 0){ 
      cat("\n")
      cat("Number of knots used: ",x$controls[8], "\n")
      cat("Value of the smoothing parameters: ", x$k0[-3], sep=" ")
      cat("\n")
    }
  }
}