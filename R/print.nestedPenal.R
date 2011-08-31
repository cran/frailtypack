
"print.nestedPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
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
	
	if (is.null(coef)){
		x$varH<-matrix(x$varH) 
		x$varHIH<-matrix(x$varHIH)
	}
 #AD:     
	if (x$typeof == 0){	 	 
		if (x$n.knots.temp < 4){
			cat("\n")
			cat("  The minimum number of knots is 4","\n")	
		}
		if (x$n.knots.temp > 20){
			cat("\n")
			cat("  The maximum number of knots is 20","\n")	
		}  
		
		if ((x$indic.Kappa2 == 0) & (x$nst == 1)){
			cat(" Kappa2 is not used  \n")
		}
	}else{
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)) cat("  The maximum number of time intervals is 20","\n")
	} 
#AD    
	if (x$istop == 1){
		if (!is.null(coef)){ 
			seH <- sqrt(diag(x$varH))[-c(1,2)]
			seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
			}	
			cat("\n")
			cat(" Nested Frailty model parameter estimates using ", "\n")

			if (x$typeof == 0){
				cat(" using a Penalized Likelihood on the hazard functions",  "\n")
			}else{
				cat(" using a Parametrical approach for the hazard functions",  "\n")	
			}
			if (x$n.strat>1) cat(" (Stratification structure used)", "\n") 
			
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

		alpha <- x$alpha
		temp <- diag(x$varH)[1]
		seH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)
		temp <- diag(x$varHIH)[1]
		seHIH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)
		
		if (x$noVar1 == 1){
			cat("\n")
			cat("    Shared Gamma Frailty model: No covariates \n")
			cat("    -------------------------- \n")
			cat("\n")
		}
		
		cat("  Frailty parameters: \n")
		if (x$typeof == 0){
			cat("   alpha  (group effect): ", alpha, " (SE(H):", 
				seH.alpha, ")", " (SE(HIH):", seHIH.alpha, ")", sep="", "\n")
		}else{
			cat("   alpha  (group effect): ", alpha, " (SE(H):", 
				seH.alpha, ")", sep="", "\n")
		}
		
		eta <- x$eta
		temp <- diag(x$varH)[2]
		seH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)
		temp <- diag(x$varHIH)[2]
		seHIH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)
		if (x$typeof == 0){
			cat("   eta (subgroup effect): ", eta, " (SE(H):", 
				seH.eta, ")", " (SE(HIH):", seHIH.eta, ")", sep="", "\n")
		}else{
			cat("   eta (subgroup effect): ", eta, " (SE(H):", 
				seH.eta, ")", sep="", "\n")
		}
		
		cat(" \n")


		if (x$typeof == 0){
			cat(paste("    penalized marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
			cat("    LCV = the approximate likelihood cross-validation criterion\n")
			cat("          in the semi parametrical case     =",x$LCV,"\n")
		}else{
			cat(paste("    marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
#			cat("    LCV = the approximate likelihood cross-validation criterion\n")
#			cat("          in the parametrical case     =",x$LCV,"\n")	
			cat("    AIC = Aikaike information Criterion     =",x$AIC,"\n")
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
		cat("    n=", x$n)
		if (length(x$na.action)){ 
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		cat("    n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter,"\n")
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)){
			cat("      Exact number of time intervals used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("    Exact number of time intervals used: ",x$nbintervR,"\n")
		 }
		if (x$typeof == 0){     
			cat("    Exact number of knots used: ", x$n.knots, "\n")
			if (x$nst == 1){
				cat("    Value of the smoothing parameter: ", x$kappa[1])
			}else{
				cat("    Value of the smoothing parameter: ", x$kappa[1], " and kappa2=", x$kappa[2], sep="")
			}
			cat(", DoF: ", formatC(-x$DoF, format="f",dig=2))
		}
	}else{
		cat("\n")
		cat(" Nested Frailty model parameter estimates using ", "\n")

		if (x$typeof == 0){
			cat(" using a Penalized Likelihood on the hazard functions",  "\n")
		}else{
			cat(" using a Parametrical approach for the hazard functions",  "\n")	
		}
		if (x$n.strat>1) cat(" (Stratification structure used)", "\n") 	
	
		if (x$noVar1 == 1){
			cat("\n")
			cat("    Nested Gamma Frailty model: No covariates \n")
			cat("    -------------------------- \n")
			cat("\n")
		}
		
		cat("\n")
		cat("    n=", x$n)
		if (length(x$na.action)){ 
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		cat("    n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter)	
	}
	cat("\n")
	invisible()
}
