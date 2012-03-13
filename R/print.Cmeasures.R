
print.Cmeasures <- function(x, digits = 3, ...) 
{
	cl <- x$call
	if (!is.null(cl)){
		cat("\n")
		cat("--------- Model ---------\n")
		dput(cl)
		cat("\n")
	}
	if(class(x)!="Cmeasures"){
		stop("The object x must be a class Cmeasures.")
	}else{
		cat("--------- Frequencies ---------\n")
		cat("\n")
		print(x$frequencies,row.names=F,digits=digits)
		cat("\n")
		if(x$Nboot > 0){
			cat("Bootstrap\n")
			 cat("     Number of bootstrap steps:",x$Nboot,"\n")
			 cat("     Number of non convergence (excluded):",x$Nbproblem,"\n")
			cat("\n")
		}
		cat("\n")
		if(x$Frailty){
			cat("--------- Concordance Probability Estimation (conditional) ---------\n")
		}else{
			cat("In Cox proportional hazards models, only marginal values\n")
			cat("of concordance probability estimation are proposed.\n")
			cat("\n")
			cat("--------- Concordance Probability Estimation ---------\n")
		}
		cat("\n")	
		print(round(x$CPEcond,digits=digits),row.names=T)
		cat("\n")
		cat("\n")
		if((x$marginal==1) & (x$Frailty)){
			cat("--------- Concordance Probability Estimation (marginal) ---------\n")	
			cat("\n")
			print(round(x$CPEmarg,digits=digits),row.names=T)
			cat("\n")
			cat("\n")
		}
		if(x$cindex==1){	
			if(x$Frailty){
				cat("--------- C-index (conditional) ---------\n")
			}else{
				cat("--------- C-index ---------\n")
			}
			cat("\n")
			print(round(x$cindexcond,digits=digits),row.names=T)	
			cat("\n")
			cat("\n")
		}
		if((x$marginal==1) & (x$cindex==1) & (x$Frailty)){
			cat("--------- C-index (marginal) ---------\n")	
			cat("\n")
			print(round(x$cindexmarg,digits=digits),row.names=T)	
		}
	}
}

