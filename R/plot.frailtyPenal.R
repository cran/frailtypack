#' Plot Method for a Shared frailty model.
#' 
#' Plots estimated baseline survival and hazard functions from an object of
#' class 'frailtyPenal'. Confidence bands are allowed.
#' 
#' 
#' @aliases plot.frailtyPenal lines.frailtyPenal
#' @usage
#' 
#' \method{plot}{frailtyPenal}(x, type.plot = "Hazard", conf.bands=TRUE,
#' pos.legend = "topright", cex.legend=0.7, main, color=2, median=TRUE, Xlab = "Time", Ylab
#' = "Hazard function", ...)
#' @param x A shared frailty model, i.e. a \code{frailtyPenal} class object
#' (output from calling \code{frailtyPenal} function).
#' @param type.plot a character string specifying the type of curve. Possible
#' value are "Hazard", or "Survival". The default is "Hazard". Only the first
#' letters are required, e.g "Haz", "Su"
#' @param conf.bands Logical value. Determines whether confidence bands will be
#' plotted.  The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#' 'par("cex")'. Default is 0.7
#' @param main title of plot
#' @param color color of the curve (integer)
#' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param ... other unused arguments
#' @return Print a plot of a shared frailty model.
#' @seealso \code{\link{frailtyPenal}}
#' @keywords file
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(readmission)
#' 
#' ###--- Shared frailty model ---###
#' 
#' modSha <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id),
#' n.knots=10,kappa=10000,data=readmission,hazard="Splines")
#' 
#' plot(modSha,type="Survival",conf=FALSE)
#' 
#' ###--- Cox proportional hazard model ---###
#' 
#' modCox <- frailtyPenal(Surv(time,event)~as.factor(dukes),n.knots=10,
#' kappa=10000,data=readmission,hazard="Splines")
#' 
#' plot(modCox)
#' 
#' #-- no confidence bands
#' plot(modSha,conf.bands=FALSE)
#' plot(modCox,conf.bands=FALSE)
#' 
#' }
#' 
#' 
"plot.frailtyPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, median=TRUE, Xlab = "Time", Ylab = "Hazard function", ...)
{
  plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be 'Hazard' or 'Survival'")
  }	
  oldpar <- par(no.readonly = TRUE) 
	on.exit(par(oldpar))
  if(missing(main)) main<-""
  if(!is.null(x$family)){
    x$family2 = x$family
  }else{
    x$family2 = 0
  }
  
  if(plot.type==1){ # hazard
    # browser()
    if(conf.bands){
      if( (is.null(x$family)) | (x$family2 %in% c(0,1,2,4)) ){
        par(mfrow=c(1,1))
        matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        for (i in (1:x$n.strat)[-1]) matlines(x$x[,i], x$lam[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
      }else{
        par(mfrow=c(1,2))
        matplot(x$x[,1], x$lam[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,
                ylab="Baseline hazard function", main=main, ...)
        trapz <- function(x,y){
          idx = 2:length(x)
          return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
        }
        nblignes = nrow(x$lam[,,1])-1
        matcumul = matrix(NA, nrow = nblignes, ncol = 3)
        abs = x$x
        ord = x$lam[,,1]
        for(j in 1:nblignes){
          matcumul[j, ] = c(
            trapz(abs[1:(j+1)], ord[1:(j+1), 1]),
            trapz(abs[1:(j+1)], ord[1:(j+1), 2]),
            trapz(abs[1:(j+1)], ord[1:(j+1), 3])
          )
        }
        matplot(x=abs[-c(1)],
                y=cbind(matcumul[,(1:3)]),
                col=color, type="l", lty=c(1,2,2),
                main=main,
                xlab=Xlab, 
                ylab="Cumulative baseline hazard function"
        )
      }
    }else{
      if( (is.null(x$family)) | (x$family2 %in% c(0,1,2,4)) ){
        par(mfrow=c(1,1))
        plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        for (i in (1:x$n.strat)[-1]) lines(x$x[,i], x$lam[,1,i], col=color+(i-1), type="l", lty=1, ...)
      }else{
        par(mfrow=c(1,2))
        plot(x$x[,1], x$lam[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab="Baseline hazard function", main=main, ...)
        trapz <- function(x,y){
          idx = 2:length(x)
          return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
        }
        nblignes = nrow(x$lam[,,1])-1
        matcumul = matrix(NA, nrow = nblignes, ncol = 3)
        abs = x$x
        ord = x$lam[,,1]
        for(j in 1:nblignes){
          matcumul[j, ] = c(
            trapz(abs[1:(j+1)], ord[1:(j+1), 1]),
            trapz(abs[1:(j+1)], ord[1:(j+1), 2]),
            trapz(abs[1:(j+1)], ord[1:(j+1), 3])
          )
        }
        plot(x=abs[-c(1)], y=matcumul[,1], col=color, type="l", lty=1, xlab=Xlab,ylab="Cumulative baseline hazard function", main=main, ...)
        # matplot(x=abs[-c(1)],
        #         y=cbind(matcumul[,(1:3)]),
        #         col=color, type="l", lty=c(1,2,2),
        #         main=main,
        #         xlab=Xlab, 
        #         ylab="Cumulative baseline hazard function"
        # )
      }
      
    }
    
  }else{ # survival
    if (missing (Ylab)) Ylab <- "Baseline survival function"
    if (x$typeof == 0){
      if (conf.bands){
        matplot(x$x[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        for (i in (1:x$n.strat)[-1]) matlines(x$x[,i], x$surv[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$x[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        for (i in (1:x$n.strat)[-1]) lines(x$x[,i], x$surv[,1,i], col=color+(i-1), type="l", lty=1, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }else{
      if (conf.bands){
        matplot(x$xSu[,1], x$surv[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main,...)
        for (i in (1:x$n.strat)[-1]) matlines(x$xSu[,i], x$surv[,,i], col=color+(i-1), type="l", lty=c(1,2,2), ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$xSu[,1], x$surv[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        for (i in (1:x$n.strat)[-1]) lines(x$xSu[,i], x$surv[,1,i], col=color+(i-1), type="l", lty=1, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }		
  }
  
  if (x$n.strat > 1) legend(pos.legend, paste("strata =",1:x$n.strat), lty=1, col=color+(1:x$n.strat-1), cex=cex.legend, ...)
  #else legend(pos.legend, c("event"), lty=1, col=color, cex=cex.legend, ...)
  
  return(invisible())
}
