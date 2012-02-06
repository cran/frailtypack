
"print.jointPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
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
	
	if (is.null(coef)){
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
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)) cat("  The maximum number of time intervals nb.int1 is 20","\n")
		if ((x$typeof == 1) & (x$indic.nb.int2 == 1)) cat("  The maximum number of time intervals nb.int2 is 20","\n")
	}     
#AD    
	if (x$istop == 1){
		if (!is.null(coef)){ 
			seH <- sqrt(diag(x$varH))[-c(1,2)]
			seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
				if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
				if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
			}
			cat("\n")
			cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}
			if (x$typeof == 0){
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test_d==1){
					dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "SE coef (HIH)", "z", "p"))
			}else{
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test_d==1){
					dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "z", "p"))
			}
	
		
			if (x$noVar1 == 0){
				cat("\n")
				cat("Recurrences:\n")
				cat("------------- \n")
				prmatrix(tmp[1:x$nvar[1], ,drop=FALSE])
				if(x$global_chisq.test==1){
					cat("\n")
					prmatrix(tmpwald)
				}
			}
			cat("\n")
		
			if (x$noVar2 == 0){	
				cat("Terminal event:\n")
				cat("---------------- \n")
				prmatrix(tmp[-c(1:x$nvar[1]), ,drop=FALSE])
				if(x$global_chisq.test_d==1){
					cat("\n")
					prmatrix(tmpwalddc)
				}
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
		if (x$typeof == 0){
			cat("   theta (variance of Frailties, Z):", theta, "(SE (H):",
				seH.theta, ")", "(SE (HIH):", seHIH.theta, ")", "\n")
			cat("   alpha (Z^alpha for terminal event):", x$alpha, "(SE (H):",
				sqrt(diag(x$varH))[2], ")", "(SE (HIH):", sqrt(diag(x$varHIH))[2],")", "\n")
			cat(" \n")
		}else{
			cat("   theta (variance of Frailties, Z):", theta, "(SE (H):",
				seH.theta, ")", ")", "\n")
			cat("   alpha (Z^alpha for terminal event):", x$alpha, "(SE (H):",
				sqrt(diag(x$varH))[2], ")", "\n")
			cat(" \n")	
		}
		if (x$typeof == 0){
			cat(paste("   penalized marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
			cat("   LCV = the approximate likelihood cross-validation criterion\n")
			cat("         in the semi parametric case     =",x$LCV,"\n")
		}else{
			cat(paste("   marginal log-likelihood =", round(x$logLik,2)))
			cat("\n")
#			cat("   LCV = the approximate likelihood cross-validation criterion\n")
#			cat("         in the parametric case     =",x$LCV,"\n")	
			cat("   AIC = Aikaike information Criterion     =",x$AIC,"\n")
			cat("\n")
			cat("The expression of the Aikaike Criterion is:","\n")
			cat("        'AIC = (1/n)[np - l(.)]'","\n")
			if (x$typeof == 2){
				cat("\n")
				cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
				cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
				cat("\n")
				cat("The expression of the Weibull hazard function is:","\n")
				cat("        'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'","\n")
				cat("The expression of the Weibull survival function is:","\n")
				cat("        'S(t) = exp[- (t/scale)^shape]'")
				cat("\n")
			}
		}
#AD:
		cat("\n")
		cat("   n=", x$n)
		if (length(x$na.action)){
			cat("      (", length(x$na.action), " observation deleted due to missing) \n")
		}else{ 
			cat("\n")
		}
		cat("   n recurrent events=", x$n.event, " n groups=", x$groups)
		cat("\n")
		cat("   n terminal events=", x$n.deaths)
		cat("\n")
		cat("   number of iterations: ", x$n.iter,"\n")
		
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)){
			cat("   Exact number of time intervals nb.int1 used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("   Exact number of time intervals nb.int1 used: ",x$nbintervR,"\n")
		 }
		if ((x$typeof == 1) & (x$indic.nb.int2 == 1)){
			cat("   Exact number of time intervals nb.int2 used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("   Exact number of time intervals nb.int2 used: ",x$nbintervDC,"\n")
		 }
		 
		if (x$typeof == 0){ 
			cat("\n")
			cat("   Exact number of knots used: ", x$n.knots, "\n")
			cat("   Value of the smoothing parameters: kappa1=", x$kappa[1], " and kappa2=", x$kappa[2], sep="")
			cat("\n")
		}
	}else{
		if (!is.null(coef)){ 
			cat("\n")
			cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}	
			
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
			
			cat("\n")
			cat("   n=", x$n)
			if (length(x$na.action)){
				cat("      (", length(x$na.action), " observation deleted due to missing) \n")
			}else{ 
				cat("\n")
			}
			cat("   n recurrent events=", x$n.event, " n groups=", x$groups)
			cat("\n")
			cat("   n terminal events=", x$n.deaths)
			cat("\n")
		}
	}
    invisible()
}
