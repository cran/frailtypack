#' Summary method for a joint competing risks midel
#' 
#' Prints a short summary of parameter estimates of a joint competing risks model
#' or more generally an object of class 'jointRecCompet'.
#' 
#' 
#' @usage \method{summary}{jointRecCompet}(object, digits = max(options()$digits - 4, 6),
#' ...)
#' @param object the result of a call to the jointRecCompet function
#' @param digits number of digits to print
#' @param \dots other unused arguments
#' @return
#' 
#' Print, separately for each type of event (Recurrent, Terminal1 and
#' Terminal2), the parameter estimates of the survival or hazard functions.
#' @seealso \code{\link{jointRecCompet}}
#' @keywords methods jointRecCompet
#' @export
"summary.jointRecCompet"<-function(object, digits = max(options()$digits - 4, 6), ...){
  x<-object
  level=.95
  z=abs(qnorm((1-level)/2))
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$summary.table$Estimate
  nvar <- length(x$coef)
  if (x$critCV[2] == 1){
    if (!is.null(coef)){
      
      seH <- sqrt(diag(x$varH.Estimate))

      
      if(x$noVarEvent[1] == 0){
        
        cat("\n")
        cat("Recurrent event:\n")
        cat("------------ \n")
        xx=grep("Recurrent", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Recurrent: ","",tmp$Parameter)
        tmp$H0=gsub("Recurrent: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          xx=grep("Weibull", tmp[,1])
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-xx,]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)
        }else{
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-xx,]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)

        }
        cat('\n')
      }
      
      
      if(x$noVarEvent[2] == 0){	
        cat("Terminal event 1:\n")
        cat("------------- \n")
        xx=grep("Terminal1", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Terminal1: ","",tmp$Parameter)
        tmp$H0=gsub("Terminal1: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          xx1=grep("Weibull", tmp[,1])
          xx2=grep("Alpha", tmp[,1])
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-c(xx1,xx2),]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)
        }else{
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          xx=grep("Alpha", tmp[,1])
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-xx,]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)
        }
        cat("\n")
      }
      
      
      if(x$noVarEvent[4]== 0){	
        cat("Terminal event 2:\n")
        cat("------------- \n")
        xx=grep("Terminal2", x$summary.table$Parameter)
        tmp = x$summary.table[xx,]
        tmp$Parameter=gsub("Terminal2: ","",tmp$Parameter)
        tmp$H0=gsub("Terminal2: ","",tmp$H0)
        tmp = tmp[,-c(2,3,6,7,9)]
        tmp$exp = exp(tmp$Estimate)
        tmp$z = tmp$Estimate/tmp$Estimate.SE
        tmp = tmp[,c(1,2,5,3,6,4)]
        
        if(x$controls[3] != 0){
          tmp$Parameter[1:2] = c('Weibull: Shape',"Weibull: Scale")
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          xx1=grep("Weibull", tmp[,1])
          xx2=grep("Alpha", tmp[,1])
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-c(xx1,xx2),]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            i=1
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)
        }else{
          colnames(tmp) = c("","coef","HR","SE coef","z","p")
          xx=grep("Alpha", tmp[,1])
          tmp$lo = exp(tmp$coef-z*tmp$`SE coef`)
          tmp$up = exp(tmp$coef+z*tmp$`SE coef`)
          tmp$coef=NULL
          tmp$`SE coef`=NULL 
          tmp$z=NULL
          tmp$p=NULL
          tmp=tmp[-xx,]
          tp=data.frame("var"='NULL',"hr"='NULL'," 95% CI"='NULL')
          for(i in 1:nrow(tmp)){
            temp=tmp[i,]
            var=as.character(temp[1])
            hr=formatC(as.numeric(tmp[2]),2,6,format="f")
            lo=formatC(as.numeric(tmp[3]),2,6,format="f")
            up=formatC(as.numeric(tmp[4]),2,6,format="f")
            ci=paste("(",lo," ; ",up," )",sep="")
            tp=rbind(tp,c(var,hr,ci))
          }
          colnames(tp) = c("",'hr',' 95%   CI    ')
          tp=tp[2:nrow(tp),]
          print(tp,,drop=FALSE,row.names = FALSE)
        }
        
        cat("\n")
      }
      
      
    }
    
    #sigma
    tmp = x$summary.table
    xsig=grep("Sigma",tmp$Parameter)
    sigma <- tmp[xsig,"Estimate"]
    sesigma <- tmp[xsig,"Estimate.SE"]
    
    
    #AD:
    if (x$noVarEvent[1] == 1){
      cat("\n")
      cat("    Recurrent event: No covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVarEvent[2] == 1){
      cat("\n")
      cat("    Terminal event 1: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    if (x$noVarEvent[4] == 1){
      cat("\n")
      cat("    Terminal event 2: No covariates \n")
      cat("    ----------- \n")
    }
  }
}