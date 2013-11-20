
"plot.predJoint" <- function (x, conf.bands=FALSE, pos.legend="topleft", cex.legend=0.7, ...)
{
	if (x$npred >10) {
		warning("For a better overview only predictions for a maximum of 10 subjects are plotted (the 10 first ones)")
		x$npred <- 10
	}
	
	if ((conf.bands) & (!x$icproba)) stop("Confidence interval was not calculated with 'prediction'")
	
	if (x$npred <= 5) par(mfrow=c(1,x$npred))
	else par(mfrow=c(2,ceiling(x$npred/2)))
	
	if (!x$intcens){ #(x$joint.clust==1){
		for (i in 1:(x$npred)) {
			if (conf.bands) {
				matplot(x$time,cbind(x$proba1[i,],x$probalow1[i,],x$probahigh1[i,]),col="blue",type="l",lty=c(1,2,2),xlab="Time",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=c(0,1))
				matlines(x$time,cbind(x$proba2[i,],x$probalow2[i,],x$probahigh2[i,]),col="red",type="l",lty=c(1,2,2))
				matlines(x$time,cbind(x$proba3[i,],x$probalow3[i,],x$probahigh3[i,]),col="green",type="l",lty=c(1,2,2))
			}else{
				plot(x$time,x$proba1[i,],col="blue",type="l",lty=c(1,2,2),xlab="Time",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=c(0,1))
				lines(x$time,x$proba2[i,],col="red",type="l",lty=c(1,2,2))
				lines(x$time,x$proba3[i,],col="green",type="l",lty=c(1,2,2))
			}
			if (i==1) legend(pos.legend, c("p1: exactly j recurrences","p2: at least j recurrences","p3: ignoring recurrences"), lty=1, col=c("blue","red","green"), cex=cex.legend)
		}
	}else{
		for (i in 1:(x$npred)) {
			if (conf.bands) {
				matplot(x$time,cbind(x$proba2[i,],x$probalow2[i,],x$probahigh2[i,]),col="red",type="l",lty=c(1,2,2),xlab="Time",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=c(0,1))
			}else{
				plot(x$time,x$proba2[i,],col="red",type="l",lty=c(1,2,2),xlab="Time",ylab="Prediction probability of event",main=paste("id patient :",x$group[i]),ylim=c(0,1))
			}
			if (i==1) legend(pos.legend, c("probability of death"), lty=1, col=c("red"), cex=cex.legend)
		}
	}
	return(invisible())
}
