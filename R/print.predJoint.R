
print.predJoint <- function(x, digits = 3, ...) 
{
# 	cl <- x$call
# 	if (!is.null(cl)){
# 		cat("\n")
# 		cat("--------- Model ---------\n")
# 		dput(cl)
# 		cat("\n")
# 	}
	if(class(x)!="predJoint"){
		stop("The object x must be a class predJoint.")
	}else{
		if (!is.null(cl <- x$call)){
			cat("Call:\n")
			dput(cl)
			cat("\n")
		}
		if (!x$intcens) {
			cat("\n")
			cat("--------- Probability 1 (exactly j recurrences) ---------\n")
			cat("\n")
			print(x$proba1,row.names=F,digits=digits)
			
			cat("\n")
			cat("--------- Probability 2 (at least j recurrences) ---------\n")
			cat("\n")
			print(x$proba2,row.names=F,digits=digits)
			
			cat("\n")
			cat("--------- Probability 3 (only parameters) ---------\n")
			cat("\n")
			print(x$proba3,row.names=F,digits=digits)
		}else{
			cat("\n")
			cat("--------- Probability  ---------\n")
			cat("\n")
			print(x$proba2,row.names=F,digits=digits)
		}
	}
}

