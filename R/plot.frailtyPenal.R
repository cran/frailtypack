"plot.frailtyPenal" <- "plot.nestedPenal" <- "plot.additivePenal" <-
function (x, type.plot="hazard", conf.bands=TRUE, ...) 
{
  
   plot.type <- charmatch(type.plot, c("hazard", "survival"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
    }

  if(plot.type==1)
   {

    if (x$n.strat==1)
      {  
        if (conf.bands)
            matplot(x$x1, x$lam, col=1, type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
        else
            plot(x$x1, x$lam[,1], col=1, type="l", lty=c(1,2,2), xlab="Time",
                ylab="Hazard function", ...)
      }
    else
      {
       if (conf.bands)
         {
           matplot(x$x1, x$lam, col=1, type="l", lty=c(1,2,2), xlab="Time",
                  ylab="Hazard function", ...)
           matlines(x$x2, x$lam2, col=2, type="l", lty=c(1,2,2),...)
         }
       else
         {
           plot(x$x1, x$lam[,1], col=1, type="l", lty=c(1,2,2), xlab="Time",
                  ylab="Hazard function", ...)
           lines(x$x2, x$lam2[,1], col=2, type="l", lty=c(1,2,2),...)
         }
  
      legend("topleft", c("strata=1", "strata=2"), lty=2, col=c(1,2))

      } 
   }        
    
  else
   {

    if (x$n.strat==1)
      {  
        if (conf.bands)
           matplot(x$x1, x$surv, col=1, type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ...)
        else        
           plot(x$x1, x$surv[,1], col=1, type="l", lty=c(1,2,2), xlab="Time",
                ylab="Baseline survival function", ...)
      }
    else
      {
        if (conf.bands)
          {
            matplot(x$x1, x$surv, col=1, type="l", lty=c(1,2,2), xlab="Time",
                   ylab="Baseline survival function", ...)
            matlines(x$x2, x$surv2, col=2, type="l", lty=c(1,2,2), ...)
          }
        else 
          {
            plot(x$x1, x$surv[,1], col=1, type="l", lty=c(1,2,2), xlab="Time",
                   ylab="Baseline survival function", ...)
            lines(x$x2, x$surv2[,1], col=2, type="l", lty=c(1,2,2), ...)
          }
        
        legend("topleft", c("strata=1", "strata=2"), lty=2, col=c(1,2))

      } 
   }        

    return(invisible())
}
