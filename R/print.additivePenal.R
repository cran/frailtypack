
"print.additivePenal" <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" frailtyPenal failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- x$coef
    nvar <- length(x$coef)
 
    if (is.null(coef))
      {
        x$varH<-matrix(x$varH) 
        x$varHIH<-matrix(x$varHIH)
      }
 #AD:     
 	 
    if (x$n.knots.temp < 4){
    	cat("\n")
        cat("  The minimum number of knots is 4","\n")	
	cat("\n")
    }  
    if (x$n.knots.temp > 20){
    	cat("\n")
         cat("  The maximum number of knots is 20","\n")	
    }        
#AD    
    if (!is.null(coef)) 
      { 

# JRG sep '09
#        if (nvar>1)
#         {    
#          seH <- sqrt(diag(x$varH))
#          seHIH <- sqrt(diag(x$varHIH))
#         }
#        if (nvar==1)
#         {    
          seH <- sqrt(x$varH)
          seHIH <- sqrt(x$varHIH)
#         }

        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
        pchisq((coef/seH)^2, 1), digits - 1))

        cat("\n")
        cat("  Additive gaussian frailty model parameter estimates ", 
            "\n")
        cat("  using a Penalized Likelihood on the hazard function", 
           "\n")

        if (x$n.strat>1)  
        cat("  (Stratification structure used)", "\n")         
	
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "SE coef (H)", "SE coef (HIH)", "z", "p"))
        cat("\n")
        prmatrix(tmp)
        cat("\n")
      } 

    if (x$correlation)
    cat("    Covariance (between the two frailty terms, \n                the intercept and the slope):", x$cov, "(SE:",x$varcov^0.5, ")", "\n")

#AD:
	if (x$noVar == 1){
		cat("\n")
		cat("    Additive gaussian frailty model: No covariates \n")
		cat("    ------------------------------- \n")
		cat("\n")
	}
#AD:
    if (x$rho!=-1)      
    cat("    Corresponding correlation between the two frailty terms :", x$rho, "\n")

    cat("    Variance for random intercept:", x$sigma2, "(SE (H):", 
            x$varSigma2[1]^.5, ")", "(SE (HIH):", x$varSigma2[2]^.5, ")", "\n")

    cat("    Variance for random slope:", x$tau2, "(SE (H):", 
            x$varTau2[1]^.5, ")", "(SE (HIH):", x$varTau2[2]^.5, ")", "\n")

    cat(" \n")
    cat(paste("    penalized marginal log-likelihood =", round(x$logLikPenal, 
        2)))
#AD:
    cat("\n")
    cat("    LCV = the approximate likelihood cross-validation criterion\n")
    cat("    in the semi parametric case     =",x$LCV,"\n")
#AD:
    cat("\n")
    cat("    n=", x$n)
    if (length(x$na.action)) 
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    else cat("\n")
    cat("    n events=", x$n.event, " n groups=", x$groups)
    cat( "\n")
    cat("    number of iterations: ", x$n.iter)
    cat("\n")
    cat("    Exact number of knots used: ", x$n.knots, "\n")

    if (!x$cross.Val)
      {
       cat("    Value of the smoothing parameter: ", x$kappa[1])
       if (x$n.strat==2)
        cat(" ", x$kappa[2])
      }

    if (x$cross.Val)
      {
       if (is.null(x$theta))
         cat("    Smoothing parameter estimated by Cross validation: ", x$kappa[1])  
       else 
         {
           cat("    Best smoothing parameter estimated by")
           cat("\n")
           cat("       an approximated Cross validation: ", x$kappa[1])
         } 
      }
    
    cat(", DoF: ", formatC(-x$DoF, format="f",dig=2))


    cat("\n")
    invisible()
}
