#' Summary of parameter estimates of a Weibull Illness-Death model with (or without) shared frailty between transitions.
#' 
#'  This function returns hazard rations (HR) and its confidence intervals
#'
#'
#' @aliases summary.frailtyIllnessDeath print.summary.frailtyIllnessDeath
#' @usage \method{summary}{frailtyIllnessDeath}(object, level = 0.95, len = 6, d = 2,
#'                                  lab="hr",...)
#' @param object output from a call to frailtyIllnessDeath.
#' @param level significance level of confidence interval. Default is 95\%.
#' @param len the total field width. Default is 6.
#' @param d the desired number of digits after the decimal point. Default of 6
#'     digits is used.
#' @param lab label of printed results.
#' @param ... other unused arguments.
#' @return
#'   Prints HR and its confidence intervals. Confidence level is allowed
#'   (level argument).


#' 
#' @examples

#'
#'   \donttest{
#'
#'
#'
#'     ###--- Semi-Markovian Weibull Illness-Death model with left truncation ---###
#'
#'     data(Paq810)
#'
#'     ModIllnessDeath_LeftTrunc <- frailtyIllnessDeath(formula = Surv(e,r,dementia) ~ gender+certif,
#'     formula.terminalEvent = Surv(t,death) ~ gender+certif  ,
#'     data=Paq810, print.info = FALSE, maxit=100)
#'
#'     #-- confidence interval at 95\% level (default)
#'
#'     summary(ModIllnessDeath_LeftTrunc)
#'
#'     #-- confidence interval at 99\% level
#'
#'     summary(ModIllnessDeath_LeftTrunc,level=0.99)
#'
#'   }


#' @seealso
#'   \code{\link{frailtyIllnessDeath}}
#' 
#' 
#' @keywords methods
#' @export
"summary.frailtyIllnessDeath" <- function(object,level=.95, len=6, d=2, lab="hr",...)
{

  
    x <- object
    if (!inherits(x, "frailtyIllnessDeath")) 
      stop("Object must be of class 'frailtyIllnessDeath'")
    
    nvar <- length(x$coef)
    
    if (is.null(x$coef)){
      if(x$Frailty==TRUE){
      cat("Weibull Illness-Death model with shared frailty: No covariates and no confidence interval\n")
      }
      
      if(x$Frailty==FALSE){
        cat("Weibull Illness-Death model: No covariates and no confidence interval\n")
      }
    }
    
    if (!is.null(x$coef)){
        
        z<-abs(qnorm((1-level)/2))
        co <- x$coef
        if (x$Frailty) {
          vcov_coef <- x$vcov[8:x$npar, 8:x$npar]
        } else {
          vcov_coef <- x$vcov[7:x$npar, 7:x$npar]
        }
        if(is.matrix(vcov_coef)){
          se <- sqrt(diag(vcov_coef))#[-1]
        }else{
          se <- sqrt(vcov_coef)
        }
        or <- exp(co)
        li <- exp(co-z * se)
        ls <- exp(co+z * se)
        r <- cbind(or, li, ls)
        dimnames(r) <- list(names(co), c(lab, paste(level*100,"%",sep=""), "C.I."))
        
        n<-r
        
        dd <- dim(n)
        n[n > 999.99] <- Inf
        a <- formatC(n, d, len,format="f")
        
        dim(a) <- dd
        
        if(length(dd) == 1){
          dd<-c(1,dd)
          dim(a)<-dd
          lab<-" "
        }else{
          lab <- dimnames(n)[[1]]
        }
        
        mx <- max(nchar(lab)) + 1
        cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
        
        for(i in (1):dd[1]) {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
          cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
        }
      }
    }

 

