#' Plot Method for a joint competing risk model with one recurrent event and 
#' two terminal events.
#' 
#' Plots of estimated baseline survival and hazard functions of joint competing 
#' recurrent model (output from an object of class 'jointRecCompet') 
#' for each type of event (recurrent and the two terminal events). 
#' Confidence intervals are allowed.
#'  
#' @usage \method{plot}{jointRecCompet}(x, event = "All", type.plot = "Hazard",
#' conf.bands = FALSE, pos.legend = "topright", cex.legend = 0.7, ylim, main,
#' color1="red", color2="blue", colorEnd="green", median=TRUE, Xlab = "Time", 
#' Ylab = "Hazard function", ...)
#' @param x A joint competing risk model, i.e. an object of class
#' \code{jointRecCompet} (output from calling \code{jointRecCompet} function).
#' @param event a character string specifying the type of outcome. Possible
#' value are "Recurrent", "Terminal1",  "Terminal2", or "All". The default is
#' "All".
#' @param type.plot a character string specifying the type of curve. Possible
#' value are "Hazard", or "Survival". The default is "Hazard". Only the first
#' words are required, e.g "Haz", "Su"
#' @param conf.bands logical value. Determines whether confidence intervals
#' will be plotted. The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#' 'par("cex")'. Default is 0.7
#' @param ylim y-axis limits
#' @param main plot title
#' @param color1 curve color for recurrent event of type 1 (integer or color
#' name in quotation marks)
#' @param color2 curve color for recurrent event of type 2 (integer or color
#' name in quotation marks)
#' @param colorEnd curve color for terminal event (integer or color name in
#' quotation marks)
#' @param median Logical value. Determines whether survival median will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Hazard function"'
#' @param \dots Other graphical parameters
#' @return Print a plot of the baseline survival or hazard functions for each
#' type of event or both with the confidence intervals or not (conf.bands
#' argument)
#' @seealso \code{\link{jointRecCompet}}
##' @export
#' @keywords methods
"plot.jointRecCompet"<-function (x, event="All", type.plot="Hazard", 
            conf.bands=FALSE, pos.legend="topright", 
            cex.legend=0.7, ylim, main, color1="red", 
            color2="blue", colorEnd="green", median=TRUE, 
            Xlab = "Time", Ylab = "Hazard function",...) 
  {
    
    event.type <- charmatch(event, c("All", "Recurrent", "Terminal1", "Terminal2"), nomatch = 0)
    if (event.type == 0) {
      stop("event must be 'All', 'Recurrent', 'Terminal1' or 'Terminal2'")
    }
    
    plot.type <- charmatch(type.plot, c("Hazard", "Survival"), 
                           nomatch = 0)
    if (plot.type == 0) {
      stop("estimator must be 'Hazard' or 'Survival'")
    }
    
    if(missing(main))
      main<-"" 
    
    x$lam = apply(x$lam,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    x$lam2 = apply(x$lam2,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    x$lam3 = apply(x$lam3,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    
    x$surv = apply(x$surv,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    x$surv2 = apply(x$surv2,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    x$surv3 = apply(x$surv3,2,function(v) sapply(v,function(x) ifelse(x==Inf,NA,x)))
    
    if (event.type==1){
      if(plot.type==1){
        if (missing(ylim)){
          yymax<-max(c(x$lam, x$lam2, x$lam3),na.rm=TRUE)
          yymin<-min(c(x$lam, x$lam2, x$lam3),na.rm=TRUE)
        }else{
          yymax<-ylim[2] 
          yymin<-ylim[1]
        }
        if (conf.bands){
          matplot(x$x1, x$lam, col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main, ...)
          matlines(x$x2, x$lam2, col=color2, type="l", lty=c(1,2,2), ...)
          matlines(x$x3, x$lam3, col=colorEnd, type="l", lty=c(1,2,2), ...)
        }else{
          plot(x$x1, x$lam[,1], col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main,...)
          lines(x$x2, x$lam2[,1], col=color2, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ...)
          lines(x$x3, x$lam3[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ...)
        } 
      }else{
        if (missing(Ylab)) Ylab <- "Baseline survival function"
        
        if (missing(ylim)){
          yymax<-1
          yymin<-0
        }else{
          yymax<-ylim[2] 
          yymin<-ylim[1]
        }
          if (conf.bands){
            matplot(x$x1, x$surv1, col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main,...)
            matlines(x$x2, x$surv2, col=color2, type="l", lty=c(1,2,2), ...)
            matlines(x$xEnd, x$survEnd, col=colorEnd, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{        
            plot(x$x1, x$surv[,1], col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=c(yymin,yymax), main=main,...)
            lines(x$x2, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), ...)
            lines(x$x3, x$surv3[,1], col=colorEnd, type="l", lty=c(1,2,2), ...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }

      }        
      legend(pos.legend, c("Recurrent event","Terminal event 1","Terminal event 2"), lty=c(1,1),col=c(color1,color2,colorEnd), xjust=1, cex=cex.legend, ...)
      
    }    
    if (event.type==2){
      
      if (missing(ylim)) ylim <- c(0,1)
      
      if(plot.type==1){
        if (conf.bands){
          matplot(x$x1, x$lam, col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        }else{
          plot(x$x1, x$lam1[,1], col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        } 
      }else{	
        if (missing(Ylab)) Ylab <- "Baseline survival function"
          if (conf.bands){
            matplot(x$x1, x$surv, col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$x1, x$surv[,1], col=color1, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
      }        
      legend(pos.legend, c("Recurrent event"), lty=c(1),col=c(color1), xjust=1, cex=cex.legend, ...)
    }
    if (event.type==3){
      
      if (missing(ylim))ylim <- c(0,1)
      
      if(plot.type==1){
        
        if (conf.bands){
          matplot(x$x1, x$lam2, col=color2, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        }else{
          plot(x$x1, x$lam2[,1], col=color2, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        } 
      }else{
        if (missing(Ylab)) Ylab <- "Baseline survival function"
          if (conf.bands){
            matplot(x$x1, x$surv2, col=color2, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{        
            plot(x$x1, x$surv2[,1], col=color2, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
      }        
      legend(pos.legend, c("Terminal event 1"), lty=c(1), col=c(color2), xjust=1, cex=cex.legend,...)
    }  
    if (event.type==4){
      
      if (missing(ylim)) ylim <- c(0,1)
      
      if(plot.type==1){
        if (conf.bands){
          matplot(x$x3, x$lam3, col=colorEnd, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        }else{
          plot(x$x3, x$lam3[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
        } 
      }else{
        if (missing(Ylab)) Ylab <- "Baseline survival function"
          if (conf.bands){
            matplot(x$x3, x$surv3, col=colorEnd, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }else{
            plot(x$x3, x$surv3[,1], col=colorEnd, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, ylim=ylim, main=main,...)
            if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
          }
      }        
      legend(pos.legend, c("Terminal event 2"), lty=c(1),col=c(colorEnd), xjust=1, cex=cex.legend, ...)
    }
    return(invisible())
  }