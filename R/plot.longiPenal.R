#' Plot Method for a joint model for longitudinal data and a terminal event.
#' 
#' Plots estimated baseline survival and hazard functions for a terminal
#' outcome from an object of class 'longiPenal'. If available, plot the estimated
#' quantities related to a mediation analysis. Confidence bands are allowed.
#' 
#' 
#' @aliases plot.longiPenal lines.longiPenal
#' @usage
#' 
#' \method{plot}{longiPenal}(x, type.plot = "Hazard",plot.mediation="All", 
#' conf.bands=TRUE,pos.legend= "topright", cex.legend=0.7, main, color, 
#' median=TRUE, Xlab = "Time", Ylab = "Hazard function", ...)
#' @param x A joint model for longitudinal outcome and a terminal event, i.e. a
#' \code{longiPenal} class object (output from calling \code{longiPenal}
#' function).
#' @param type.plot a character string specifying the type of curve for the
#' terminal event. Possible value are "Hazard", or "Survival". The default is
#' "Hazard". Only the first words are required, e.g "Haz", "Su"
#' @param plot.mediation A character string specifying the desired plot.
#'  Possible values are "All", "PTE" or "Effects". The default is
#' "All" which displays both plots.
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
#' @return Print a plot for the terminal event of the joint model for a
#' longitudinal and survival data.
#' @seealso \code{\link{longiPenal}}
#' @keywords file
##' @export
#' @examples
#' 
#' 
#' \dontrun{
#' ###--- Joint model for longitudinal data and a terminal event ---###
#' 
#' data(colorectal)
#' data(colorectalLongi)
#' 
#' # Survival data preparation - only terminal events 
#' colorectalSurv <- subset(colorectal, new.lesions == 0)
#' 
#' # Baseline hazard function approximated with splines
#' # Random effects as the link function
#' 
#' model.spli.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS 
#' + prev.resection, tumor.size ~  year * treatment + age + who.PS ,
#' colorectalSurv,	data.Longi = colorectalLongi, random = c("1", "year"),
#' id = "id", link = "Random-effects", left.censoring = -3.33, 
#' n.knots = 7, kappa = 2)
#' pdf(file = "/home/agareb1/etudiants/al10/newpack/test/plot_longi.pdf")
#' 
#' # Plot the estimated baseline hazard function with the confidence intervals
#' plot(model.spli.RE)	
#' 
#' # Plot the estimated baseline hazard function with the confidence intervals
#' plot(model.spli.RE, type = "Survival")	
#' }
#' 
#' 
"plot.longiPenal" <- function (x, type.plot="Hazard", plot.mediation='All',
                               conf.bands=TRUE, pos.legend="topright", cex.legend=0.7,
                               main, color=2, median=TRUE, Xlab = "Time", 
                               Ylab = "Hazard function",...)
{
  
  plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
  plot.type.med <- charmatch(plot.mediation, c("All", "PTE",'Effects'),nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be Hazard or Survival")
  }
  
  if(plot.type.med == 0){
    stop(" 'plot.mediation' must be either 'All', 'PTE' or 'Effects'.")
  }
  
  if(missing(main))
    main<-""
  
  if(plot.type==1){ # hazard
    
    if(conf.bands){
      matplot(x$xD[-1,1], x$lamD[-1,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
    }else{
      plot(x$xD[-1,1], x$lamD[-1,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
    }
    
  }else{ # survival
    if (missing(Ylab)) Ylab <- "Baseline survival function"
    if (x$typeof == 0){
      if (conf.bands){
        matplot(x$xD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$xD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }else{
      if (conf.bands){
        matplot(x$xSuD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }else{
        plot(x$xSuD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=main, ...)
        if (median){abline(a=0.5,b=0,cex=.5,col=1,lty=3)}
      }
    }
    
  }
  
  legend(pos.legend, c("event"), lty=1, col=color, cex=cex.legend, ...)
  
  
  if(!is.null(x$mediation)){
    if(plot.mediation=="PTE"){
      x<-x$mediation
      if(length(x)==5 & conf.bands){
        data.pte<-x$data.pte
        pte.ci=x$PTE.ci
        nie.ci=x$NIE.ci
        nde.ci = x$NDE.ci
        te.ci =x$TE.ci
        #plot r(t)
        ymin<-ifelse(min(pte.ci$lower,na.rm=T)<0,
                     1.2*min(pte.ci$lower,na.rm=T),
                     0.8*min(pte.ci$yuppower,na.rm=T))
        ymax<-ifelse(max(pte.ci$upper,na.rm = T)<0,
                     0.8*max(pte.ci$upper,na.rm=T),
                     1.2*max(pte.ci$upper,na.rm=T))
        invisible(readline(prompt="Press [Enter] to continue:"))
        plot(x<-data.pte$Time,y<-data.pte$PTE,type='l',col="black",
             xlab="Time",ylab='PTE',main="Estimated PTE with its 95% confidence band",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-pte.ci$upper,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-pte.ci$lower,type='l',lty=2,col="black",...)
      }else{
        #plot without confidence bands
        data.pte<-x$data.pte
        #plot r(t)
        invisible(readline(prompt="Press [Enter] to continue:"))
        plot(x<-data.pte$Time,y<-data.pte$PTE,type='l',col="black",
             xlab="Time",ylab='PTE',main="Estimated PTE",...)
      }
    }
    if(plot.mediation=="Effects"){
      x<-x$mediation
      if(length(x)==5 & conf.bands){
        data.pte<-x$data.pte
        pte.ci=x$PTE.ci
        nie.ci=x$NIE.ci
        nde.ci = x$NDE.ci
        te.ci =x$TE.ci
        #plot effects
        invisible(readline(prompt="Press [Enter] to continue:"))
        miny <- min(te.ci$lower,nie.ci$lower,nde.ci$lower,na.rm=T)
        maxy <-max(te.ci$upper,nie.ci$upper,nde.ci$upper,na.rm=T)
        ymin<-ifelse(miny<0,1.2*miny,0.8*miny)
        ymax<-ifelse(maxy<0,0.8*maxy,1.2*maxy)
        
        plot(x<-data.pte$Time,y<-data.pte$TE,type='l',col="black",
             xlab="Time",ylab='Effects',main="Estimated natural effects with their 95% confidence bands",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-data.pte$NDE,type='l',col="green",...)
        lines(x<-data.pte$Time,y<-data.pte$NIE,type='l',col="red",...)
        
        lines(x<-data.pte$Time,y<-te.ci$lower,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-nde.ci$lower,type='l',lty=2,col="green",...)
        lines(x<-data.pte$Time,y<-nie.ci$lower,type='l',lty=2,col="red",...)
        
        lines(x<-data.pte$Time,y<-te.ci$upper,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-nde.ci$upper,type='l',lty=2,col="green",...)
        lines(x<-data.pte$Time,y<-nie.ci$upper,type='l',lty=2,col="red",...)
        
        legend(pos.legend,legend=c("Total effect","Direct effect","Indirect effect"),
               col=c("black","green","red"),lty=1)
      }else{
        #plot without confidence bands
        data.pte<-x$data.pte
        #plot effect
        invisible(readline(prompt="Press [Enter] to continue:"))
        
        ymin<-ifelse(min(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     1.2*min(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     0.8*min(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        ymax<-ifelse(max(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     0.8*max(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     1.2*max(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        plot(x<-data.pte$Time,y<-data.pte$TE,type='l',col="black",
             xlab="Time",ylab='Effects',main="Estimated natural effects",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-data.pte$NIE,type='l',col="green",...)
        lines(x<-data.pte$Time,y<-data.pte$NDE,type='l',col="red",...)
        legend(pos.legend,legend=c("Total effect","Direct effect","Indirect effect"),
               col=c("black","green","red"),lty=1)
        
      }
    }
    if(plot.mediation=="All"){
      x<-x$mediation
      if(length(x)==5 & conf.bands){
        data.pte<-x$data.pte
        pte.ci=x$PTE.ci
        nie.ci=x$NIE.ci
        nde.ci = x$NDE.ci
        te.ci =x$TE.ci
        #plot r(t)
        ymin<-ifelse(min(pte.ci$lower,na.rm=T)<0,
                     1.2*min(pte.ci$lower,na.rm=T),
                     0.8*min(pte.ci$lower,na.rm=T))
        ymax<-ifelse(max(pte.ci$upper,na.rm = T)<0,
                     0.8*max(pte.ci$upper,na.rm=T),
                     1.2*max(pte.ci$upper,na.rm=T))
        invisible(readline(prompt="Press [Enter] to continue:"))
        plot(x<-data.pte$Time,y<-data.pte$PTE,type='l',col="black",
             xlab="Time",ylab='PTE',main="Estimated PTE with its 95% confidence band",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-pte.ci$upper,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-pte.ci$lower,type='l',lty=2,col="black",...)
        
        #plot effects
        invisible(readline(prompt="Press [Enter] to continue:"))
        miny <- min(te.ci$lower,nie.ci$lower,nde.ci$lower,na.rm=T)
        maxy <-max(te.ci$upper,nie.ci$upper,nde.ci$upper,na.rm=T)
        ymin<-ifelse(miny<0,1.2*miny,0.8*miny)
        ymax<-ifelse(maxy<0,0.8*maxy,1.2*maxy)
        
        plot(x<-data.pte$Time,y<-data.pte$TE,type='l',col="black",
             xlab="Time",ylab='Effects',main="Estimated natural effects with their 95% confidence bands",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-data.pte$NDE,type='l',col="green",...)
        lines(x<-data.pte$Time,y<-data.pte$NIE,type='l',col="red",...)
        
        lines(x<-data.pte$Time,y<-te.ci$lower,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-nde.ci$lower,type='l',lty=2,col="green",...)
        lines(x<-data.pte$Time,y<-nie.ci$lower,type='l',lty=2,col="red",...)
        
        lines(x<-data.pte$Time,y<-te.ci$upper,type='l',lty=2,col="black",...)
        lines(x<-data.pte$Time,y<-nde.ci$upper,type='l',lty=2,col="green",...)
        lines(x<-data.pte$Time,y<-nie.ci$upper,type='l',lty=2,col="red",...)
        
        legend(pos.legend,legend=c("Total effect","Direct effect","Indirect effect"),
               col=c("black","green","red"),lty=1)
      }else if(conf.bands & length(x) != 5){
        data.pte<-x$data.pte
        invisible(readline(prompt="Press [Enter] to continue:"))
        plot(x<-data.pte$Time,y<-data.pte$PTE,type='l',col="black",
             xlab="Time",ylab='PTE',main="Estimated PTE",...)
        
        
        #plot effect
        invisible(readline(prompt="Press [Enter] to continue:"))
        
        ymin<-ifelse(min(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     1.2*min(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     0.8*min(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        ymax<-ifelse(max(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     0.8*max(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     1.2*max(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        plot(x<-data.pte$Time,y<-data.pte$TE,type='l',col="black",
             xlab="Time",ylab='Effects',main="Estimated natural effects",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-data.pte$NIE,type='l',col="green",...)
        lines(x<-data.pte$Time,y<-data.pte$NDE,type='l',col="red",...)
        legend(pos.legend,legend=c("Total effect","Direct effect","Indirect effect"),
               col=c("black","green","red"),lty=1)
        
        
      }else{
        #plot without confidence bands
        data.pte<-x$data.pte
        #plot r(t)
        invisible(readline(prompt="Press [Enter] to continue:"))
        plot(x<-data.pte$Time,y<-data.pte$PTE,type='l',col="black",
             xlab="Time",ylab='PTE',main="Estimated PTE",...)
        #plot effect
        invisible(readline(prompt="Press [Enter] to continue:"))
        
        ymin<-ifelse(min(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     1.2*min(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     0.8*min(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        ymax<-ifelse(max(data.pte$TE,data.pte$NIE,data.pte$NDE)<0,
                     0.8*max(data.pte$TE,data.pte$NIE,data.pte$NDE),
                     1.2*max(data.pte$TE,data.pte$NIE,data.pte$NDE))
        
        plot(x<-data.pte$Time,y<-data.pte$TE,type='l',col="black",
             xlab="Time",ylab='Effects',main="Estimated natural effects",
             ylim=c(ymin,ymax),...)
        lines(x<-data.pte$Time,y<-data.pte$NIE,type='l',col="green",...)
        lines(x<-data.pte$Time,y<-data.pte$NDE,type='l',col="red",...)
        legend(pos.legend,legend=c("Total effect","Direct effect","Indirect effect"),
               col=c("black","green","red"),lty=1)
        
      }
    }
  }
  return(invisible())
}