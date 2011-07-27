
"print.jointPenal" <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        if (x$type == "counting") {
            cat("\n      left truncated structure used")
        }
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
      	
        seH <- sqrt(diag(x$varH))[-c(1,2)]
        seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
        pchisq((coef/seH)^2, 1), digits - 1))
        cat("\n")
        cat("  Joint gamma frailty model for recurrent and a terminal event processes", 
            "\n")
        cat("  using a Penalized Likelihood on the hazard function", 
            "\n")
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "SE coef (H)", "SE coef (HIH)", "z", "p"))

	
	if (x$noVar1 == 0){
		cat("\n")
		cat("Recurrences:\n")
		cat("------------- \n")
		prmatrix(tmp[1:x$nvar[1], ,drop=FALSE])
	}
        cat("\n")

	if (x$noVar2 == 0){	
		cat("Terminal event:\n")
		cat("---------------- \n")
		prmatrix(tmp[-c(1:x$nvar[1]), ,drop=FALSE])
		cat("\n")
	}
    }
    theta <- x$theta
    temp <- diag(x$varH)[1]
    seH.theta <- sqrt(((2 * (theta^0.5))^2) * temp)
    temp <- diag(x$varHIH)[1]
    seHIH.theta <- sqrt(((2 * (theta^0.5))^2) * temp)
#AD:
	if (x$noVar1 == 1){
		cat("\n")
		cat("    Recurrences: No covariates \n")
		cat("    ----------- \n")
	}
	
	if (x$noVar2 == 1){
		cat("\n")
		cat("    Terminal event: No covariates \n")
		cat("    -------------- \n")
		cat("\n")
	}
#AD:  
    cat(" Frailty parameters: \n")
    cat("   theta (variance of Frailties, Z):", theta, "(SE (H):",
        seH.theta, ")", "(SE (HIH):", seHIH.theta, ")", "\n")
    cat("   alpha (Z^alpha for terminal event):", x$alpha, "(SE (H):",
        sqrt(diag(x$varH))[2], ")", "(SE (HIH):", sqrt(diag(x$varHIH))[2],
        ")", "\n")
    cat(" \n")
    cat(paste("   penalized marginal log-likelihood =", round(x$logLikPenal, 
        2)))
#AD:
    cat("\n")
    cat("   LCV = the approximate likelihood cross-validation criterion\n")
    cat("   in the semi parametric case     =",x$LCV,"\n")
#AD:
    cat("\n")
    cat("   n=", x$n)
    if (length(x$na.action))
        cat("      (", length(x$na.action), " observation deleted due to missing) \n")
    else cat("\n")
    cat("   n recurrent events=", x$n.event, " n groups=", x$groups)
    cat("\n")
    cat("   n terminal events=", x$n.deaths)
    cat("\n")
    cat("   number of iterations: ", x$n.iter)
    cat("\n")
    cat("   Exact number of knots used: ", x$n.knots, "\n")
    cat("   Value of the smoothing parameters: kappa1=", x$kappa[1], " and kappa2=", x$kappa[2], sep="")
    cat("\n")
    invisible()
}
