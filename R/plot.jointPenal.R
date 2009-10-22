"plot.jointPenal" <-
function (x, event="both", type.plot="hazard", conf.bands=FALSE, ...) 
{
  
   event.type <- charmatch(event, c("both", "recurrent", "terminal"), nomatch = 0)
    if (event.type == 0) {
        stop("event must be 'both', 'recurrent' or 'terminal'")
    }


   plot.type <- charmatch(type.plot, c("hazard", "survival"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be 'hazard' or 'survival'")
    }


  if (event.type==1)
    {
     if(plot.type==1)
      {
        yymax<-max(c(x$lam, x$lam2),na.rm=TRUE)
        yymin<-min(c(x$lam, x$lam2),na.rm=TRUE)
        if (conf.bands)
          {
            matplot(x$x1, x$lam, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ylim=c(yymin,yymax), ...)
            matlines(x$x2, x$lam2, col="blue", type="l", lty=c(1,2,2), ...)
           }

        else
          {
            plot(x$x1, x$lam[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ylim=c(yymin,yymax), ...)

            lines(x$x2, x$lam2[,1], col="blue", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
          } 
      }        
    
     else
      {
        if (conf.bands)
         {
           matplot(x$x1, x$surv, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
           matlines(x$x2, x$surv2, col="blue", type="l", lty=c(1,2,2), ...)
         } 
        
        else
         {        
           plot(x$x1, x$surv[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
           lines(x$x2, x$surv2[,1], col="blue", type="l", lty=c(1,2,2), ...)
         }

      }        
        legend("topright",c("recurrent events", "terminal event"), lty=c(1,1),col=c("red","blue"),xjust=1)

   }


  if (event.type==2)
    {
     if(plot.type==1)
      {

        if (conf.bands)
          {
            matplot(x$x1, x$lam, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
           }

        else
          {
            plot(x$x1, x$lam[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
          } 
      }        
    
     else
      {
        if (conf.bands)
         {
           matplot(x$x1, x$surv, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
         } 
        
        else
         {        
           plot(x$x1, x$surv[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
         }

      }        
        legend("topright",c("recurrent events"), lty=c(1),col=c("red"),xjust=1)
   }




  if (event.type==3)
    {
     if(plot.type==1)
      {

        if (conf.bands)
          {
            matplot(x$x1, x$lam2, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
           }

        else
          {
            plot(x$x1, x$lam2[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
          } 
      }        
    
     else
      {
        if (conf.bands)
         {
           matplot(x$x1, x$surv2, col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
         } 
        
        else
         {        
           plot(x$x1, x$surv2[,1], col="red", type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ylim=c(0,1), ...)
         }

      }        
        legend("topright",c("terminal event"), lty=c(1),col=c("red"),xjust=1)
   }




    return(invisible())
}
