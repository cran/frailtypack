
print.epoce <- function(x, digits = 3, ...) 
{
	if(class(x)!="epoce"){
		stop("The object x must be a class epoce.")
	}else{
# 		if (!is.null(cl <- x$call)){
# 			cat("Call:\n")
# 			dput(cl)
# 			cat("\n")
# 		}
		if (x$new.data){
			out <- matrix(x$mpl,nrow=1,ncol=length(x$pred.times),byrow=TRUE)
			rownames(out) <- c("mpl")
		}else{
			out <- matrix(c(x$mpl,x$cvpl),nrow=2,ncol=length(x$pred.times),byrow=TRUE)
			rownames(out) <- c("mpl","cvpl")
		}
		colnames(out) <- x$pred.times
		print(out)
	}
}

