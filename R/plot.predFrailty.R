
"plot.predFrailty" <- function (x, conf.bands=FALSE, pos.legend="topleft", cex.legend=0.7, ...)
{
# 	if (x$npred >10) {
# 		warning("For a better overview only predictions for a maximum of 10 subjects are plotted (the 10 first ones)")
# 		x$npred <- 10
# 	}
# 	
# 	if ((conf.bands) & (!x$icproba)) stop("Confidence interval was not calculated with 'predictFP'")
# 	
# 	if (x$npred <= 5) par(mfrow=c(1,x$npred))
# 	else par(mfrow=c(2,ceiling(x$npred/2)))
# 	
# 	for (i in 1:(x$npred)) {
# 		if (conf.bands) {
# 			matplot(x$time,cbind(x$proba1[i,],x$probalow1[i,],x$probahigh1[i,]),col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="Probability of event",main=paste("individu",i),ylim=c(0,1))
# 			matlines(x$time,cbind(x$proba2[i,],x$probalow2[i,],x$probahigh2[i,]),col="red",type="l",lty=c(1,2,2))
# 			matlines(x$time,cbind(x$proba3[i,],x$probalow3[i,],x$probahigh3[i,]),col="green",type="l",lty=c(1,2,2))
# 		}else{
# 			plot(x$time,x$proba1[i,],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="Probability of event",main=paste("individu",i),ylim=c(0,1))
# 			lines(x$time,x$proba2[i,],col="red",type="l",lty=c(1,2,2))
# 			lines(x$time,x$proba3[i,],col="green",type="l",lty=c(1,2,2))
# 		}
# 		legend(pos.legend, c("proba1","proba2","proba3"), lty=1, col=c("blue","red","green"), cex=cex.legend, ...)
# 	}

	# par(mfrow=c(2,1))
	if (is.null(x$type)){
		title <- paste("Probability")
	}else{
		if (x$type=="marginal") title <- paste("Marginal probability")
		else title <- paste("Conditional probability")
	}
	
	if (conf.bands){
		matplot(x$time,x$pred,type="l",lty=1,xlab="Time",ylab="Prediction probability of event",main=title,ylim=c(0,1))
		matlines(x$time,x$predLow,type="l",lty=2)
		matlines(x$time,x$predHigh,type="l",lty=2)
	}else{
		matplot(x$time,x$pred,type="l",lty=1,xlab="Time",ylab="Prediction probability of event",main=title,ylim=c(0,1))
	}
	
	legend(pos.legend, paste("profile",(1:ncol(x$pred))),lty=1,bty="n",col=(1:ncol(x$pred)))

	return(invisible())
}
