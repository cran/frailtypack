
"print.frailtyPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
		if (x$type == "counting"){
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
	if (x$typeof == 0){	 
		if (x$n.knots.temp < 4){
			cat("\n")
			cat("  The minimum number of knots is 4","\n")	
			cat("\n")
		} 
		if (x$n.knots.temp > 20){
			cat("\n")
			cat("  The maximum number of knots is 20","\n")	
		}  
	}else{
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)) cat("  The maximum number of time intervals is 20","\n")
	}   
#AD     
	if (x$istop == 1){
		if (!is.null(coef)){ 
			seH <- sqrt(diag(x$varH))[-1]
			seHIH <- sqrt(diag(x$varHIH))[-1]
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
			}
		
			cat("\n")
			if (!is.null(x$theta)){
			
				cat("  Shared Gamma Frailty model parameter estimates ","\n")
				if (x$typeof == 0){
					cat("  using a Penalized Likelihood on the hazard function","\n")
				}else{
					cat("  using a Parametrical approach for the hazard function","\n")
				}
		
				if (x$n.strat>1) cat("  (Stratification structure used)", "\n")         
			}else{
				if (x$typeof == 0){
					cat("  Cox proportional hazards model parameter estimates ","\n")
					cat("  using a Penalized Likelihood on the hazard function","\n")
				}else{
					cat("  Cox proportional hazards model parameter estimates ","\n")
					cat("  using a Parametrical approach for the hazard function","\n")
				}	
				
				if (x$n.strat>1) cat("  (Stratification structure used)", "\n")         
			}

			if (x$typeof == 0){
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "SE coef (HIH)", "z", "p"))
			}else{
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "z", "p"))
			}
			cat("\n")
			prmatrix(tmp)

			cat("\n")
		} 
	
		if (!is.null(x$theta)) {
			tetha <- x$theta
			temp <- diag(x$varH)[1]
			seH <- sqrt(((2 * (tetha^0.5))^2) * temp)
			temp <- diag(x$varHIH)[1]
			seHIH <- sqrt(((2 * (tetha^0.5))^2) * temp)
#AD:
			if (x$noVar1 == 1){
				cat("\n")
				cat("    Shared Gamma Frailty model: No covariates \n")
				cat("    -------------------------- \n")
				cat("\n")
			}
#AD:		  
			if (x$typeof == 0){
				cat("    Frailty parameter, Theta:", tetha, "(SE (H):", 
				seH, ")", "(SE (HIH):", seHIH, ")", "\n")
			}else{
				cat("    Frailty parameter, Theta:", tetha, "(SE (H):", 
				seH, ")", "\n")
			}
			
        	}
		cat(" \n")     
#AD:	
		if (x$typeof == 0){
			cat(paste("      penalized marginal log-likelihood =", round(x$logLikPenal,2)))	
			cat("\n")
			cat("      LCV = the approximate likelihood cross-validation criterion\n")
			cat("            in the semi parametrical case     =",x$LCV,"\n")
		}else{
			cat(paste("      marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
#			cat("      LCV = the approximate likelihood cross-validation criterion\n")
#			cat("            in the parametrical case     =",x$LCV,"\n")	
			cat("      AIC = Aikaike information Criterion     =",x$AIC,"\n")
			cat("\n")
			cat("The expression of the Aikaike Criterion is:","\n")
			cat("        'AIC = (1/n)[np - l(.)]'","\n")
			if (x$typeof == 2){
				cat("\n")
				if (x$n.strat == 1){
					cat("      Scale for the weibull hazard function is :",round(x$shape.weib[1],2),"\n")	
					cat("      Shape for the weibull hazard function is :",round(x$scale.weib[1],2),"\n")
				}else{
					cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
					cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
				}
					cat("\n")
					cat("The expression of the Weibull hazard function is:","\n")
					cat("        'lambda(t) = shape(t^(shape-1)/(scale^shape)'","\n")
					cat("The expression of the Weibull survival function is:","\n")
					cat("        'S(t) = exp[- (t/scale)^shape]'")
					cat("\n")
			}
		}
#AD:
		cat("\n")
		cat("      n=", x$n)
		
		if (length(x$na.action)){
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		
		cat("      n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("      number of iterations: ", x$n.iter,"\n")
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)){
			cat("      Exact number of time intervals used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
		 }	
		if (x$typeof == 0){ 
			cat("\n")
			cat("      Exact number of knots used: ", x$n.knots, "\n")
	
			if (!x$cross.Val){
				cat("      Value of the smoothing parameter: ", x$kappa[1])
				if (x$n.strat==2) cat(" ", x$kappa[2])
			}
		
			if (x$cross.Val){
				if (is.null(x$theta)){
					cat("      Smoothing parameter estimated by Cross validation: ", x$kappa[1])  
				}else{ 
					cat("      Best smoothing parameter estimated by")
					cat("\n")
					cat("      an approximated Cross validation: ", x$kappa[1])
				}
			}
			cat(", DoF: ", formatC(-x$DoF, format="f",dig=2))
		}
	}else{
		if (!is.null(x$theta)){
			cat("  Shared Gamma Frailty model parameter estimates ","\n")
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}

			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")   
			if (x$noVar1 == 1){
				cat("\n")
				cat("    Shared Gamma Frailty model: No covariates \n")
				cat("    -------------------------- \n")
				cat("\n")
			}
				
		}else{
			if (x$typeof == 0){
				cat("  Cox proportional hazards model parameter estimates ","\n")
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  Cox proportional hazards model parameter estimates ","\n")
				cat("  using a Parametrical approach for the hazard function","\n")
			}	
			
			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")         
		}
		cat("\n")
		cat("      n=", x$n)
		
		if (length(x$na.action)){
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		
		cat("      n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("      number of iterations: ", x$n.iter)
	}
	cat("\n")
	invisible()
}

