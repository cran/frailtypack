#' summary of parameter estimates of a Weibull competing risks model with (or without) shared frailty between transitions.
#' 
#'  This function returns hazard rations (HR) and its confidence intervals
#'
#'
#' @aliases summary.frailtyCmprsk print.summary.frailtyCmprsk
#' @usage \method{summary}{frailtyCmprsk}(object, level = 0.95, len = 6, d = 2,
#'                                  lab="hr",...)
#' @param object output from a call to frailtyCmprsk.
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
#'     ###--- Simple Weibull competing risks model ---###
#'
#'     ###--- Weibull competing risks model with shared frailty between transitions ---###
#'
#'     data(CPRSKbmtcrr)
#'
#'     ##--- Simple Weibull competing risks model with left truncation ---##
#'
#'     modCmprskFrailty <- frailtyCmprsk(
#'       formulas = list(
#'         Surv(Age, observed_time, Status, type = "mstate") ~ Sex,
#'         ~ Sex
#'       ),
#'       data = CPRSKbmtcrr,
#'       print.info = FALSE,
#'       maxit = 100
#'     )
#'
#'     #-- Confidence interval at 95% level (default)
#'
#'     summary(modCmprskFrailty)
#'
#'     #-- Confidence interval at 99% level
#'
#'     summary(modCmprskFrailty, level = 0.99) }

#' @seealso
#'   \code{\link{frailtyCmprsk}}
#' 
#' 
#' @keywords methods
#' @export
"summary.frailtyCmprsk" <- function(object,level=.95, len=6, d=2, lab="hr",...)
{
  
  x <- object
  n_competing_events_from_formulas <- length(x$formulas)
  deb_reg <- 2*n_competing_events_from_formulas + (if (x$Frailty == TRUE) 1 else 0) + 1 # Start index for betas in x$b
  
  if (!inherits(x, "frailtyCmprsk")) 
    stop("Object must be of class 'frailtyCmprsk'")
  
  nvar <- length(x$coef)
  
  if (is.null(x$coef)){
    if(x$Frailty==TRUE){
      cat("Weibull competing risks model with shared frailty: No covariates and no confidence interval\n")
    }
    
    if(x$Frailty==FALSE){
      cat("Weibull competing risks model: No covariates and no confidence interval\n")
    }
  }
  
  if (!is.null(x$coef)){
    
    z<-abs(qnorm((1-level)/2))
    co <- x$coef
    vcov_coef <- x$vcov[deb_reg:x$npar, deb_reg:x$npar]
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



