
"plot.epoce" <- function (x, pos.legend="topright", cex.legend=0.7, ...)
{

	if (x$new.data){
		plot(x$pred.times,x$mpl,type="b",pch="X",col="blue",xlab="time",ylab="epoce")
		legend(pos.legend,c("mpl"),lty=1,col="blue")
	}else{
		plot(x$pred.times,x$cvpl,type="b",pch="X",col="red",xlab="time",ylab="epoce")
		legend(pos.legend,c("cvpl"),lty=1,col="red")
	}

	return(invisible())
}
