
"frailtyPenal" <-
function (formula, data, Frailty = TRUE, recurrentAG=FALSE, n.knots, kappa1 ,
           kappa2, maxit=350)
 {

   
    if (Frailty) 
      {
        frailty <- 1
      }

    if (!Frailty) 
      {
        frailty <- 0
      }

   
    if (missing(n.knots))
         stop("number of knots are required")       

    if (n.knots<4) n.knots<-4
    if (n.knots>20) n.knots<-20

    if (missing(kappa1))
         stop("smoothing parameter (kappa1) is required")       


    call <- match.call()
    m <- match.call(expand = FALSE)
    m$Frailty <- m$n.knots <- m$recurrentAG <- m$kappa1 <- m$kappa2 <-  m$... <- NULL
    special <- c("strata", "cluster")
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)   
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
  
    if (NROW(m) == 0) 
        stop("No (non-missing) observations")
   
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    ll <- attr(Terms, "term.labels")
    mt <- attr(m, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, m, contrasts)
    
    strats <- attr(Terms, "specials")$strata
    cluster <- attr(Terms, "specials")$cluster

   
    dropx <- NULL

    if (length(cluster)) {
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) 
            stop("Cluster can not be used in an interaction")
        cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
        dropx <- tempc$terms
        uni.cluster<-unique(cluster)
    }
    else
     {
      stop("grouping variable is needed")   
     }

    if(length(uni.cluster)==1) 
     {
      stop("grouping variable must have more than 1 level")   
     }
 
    if(length(uni.cluster)>1000) 
     {
      stop("grouping variable must have less than 1000 groups
             \n please contact to the mantainer")   
     }
    

    if (length(strats)) {
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
        uni.strat<-length(unique(strats))
      
        if (missing(kappa1))
         stop("smoothing parameter (kappa2) is required") 
    
        if (uni.strat!=2) 
          {
             stop("maximum number of strata is 2")
          }

    }
    else
      {
        uni.strat<-1
        strats <- rep(1,nrow(data))
        kappa2<-0 
      }
    
    type <- attr(Y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("Cox model doesn't support \"", type, "\" survival data", 
            sep = ""))    

    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 
        1)
    if (ncol(X) == 1) 
      {
         X<-X-1
         noVar<-1 
      }
    else
      {
         X <- X[, -1, drop = FALSE]
         noVar<-0
      }     

    nvar<-ncol(X) 

    if(nvar>15)
       stop("maximum number of variables allowed are 15. 
             \n please contact to the mantainer")

    var<-matrix(c(X),nrow=nrow(X),ncol=nvar)

    n<-nrow(X)    
    if(n>15000) 
     {
      stop("number of observations must be less than 15000 
             \n please contact to the mantainer")   
     }

    if (type=="right")
      {
        tt0 <- rep(0,n)
        tt1 <- Y[,1]
        cens <- Y[,2]
      }
    else
      {
        tt0 <- Y[,1]
        tt1 <- Y[,2]
        cens <- Y[,3]
      }

    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0

    AG<-ifelse(recurrentAG,1,0)

    ans <- .Fortran("frailpenal",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(cens.data),
                as.integer(uni.strat),
                as.integer(frailty),
                as.integer(n.knots),
                as.double(kappa1),
                as.double(kappa2),
                as.double(tt0),
                as.double(tt1),
                as.integer(cens),
                as.integer(cluster),
                as.integer(nvar),
                as.double(strats),
                as.double(var),
                as.integer(AG),
                as.integer(noVar), 
                as.integer(maxit),
                as.integer(0),
                as.double(rep(0,50)),
                as.double(matrix(0,nrow=50,ncol=50)),
                as.double(matrix(0,nrow=50,ncol=50)),
                as.double(0),
                as.double(rep(0,99)),
                as.double(matrix(0,nrow=99,ncol=3)),
                as.double(matrix(0,nrow=99,ncol=3)),
                as.double(rep(0,99)),
                as.double(matrix(0,nrow=99,ncol=3)),
                as.double(matrix(0,nrow=99,ncol=3)),
                as.integer(0),
                as.integer(0),
                as.integer(0), PACKAGE = "frailtypack"
     )    
    
    if (noVar==1) nvar<-0

    np <- ans[[19]]
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans[[31]]
    fit$logVerComPenal <- ans[[23]]
    if (Frailty) {
        fit$theta <- (ans[[20]][np - nvar])^2
    }
    if (!Frailty) {
        fit$theta <- NULL
    }
    if (noVar==1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans[[20]][(np - nvar + 1):np]
       names(fit$coef) <- colnames(X)
     }
    

    temp1 <- matrix(ans[[21]], nrow = 50, ncol = 50)[1:np, 1:np]
    temp2 <- matrix(ans[[22]], nrow = 50, ncol = 50)[1:np, 1:np]
    fit$varH <- temp1[(np - nvar):np, (np - nvar):np]
    fit$varHIH <- temp2[(np - nvar):np, (np - nvar):np]
    fit$formula <- formula(Terms)
    fit$x1 <- ans[[24]]
    fit$lam <- matrix(ans[[25]], nrow = 99, ncol = 3)
    fit$surv <- matrix(ans[[26]], nrow = 99, ncol = 3)
    fit$x2 <- ans[[27]]
    fit$lam2 <- matrix(ans[[28]], nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans[[29]], nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots 
    fit$n.iter <- ans[[30]]
    
    if(ans[[32]]==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    class(fit) <- "frailtyPenal"
    fit

    
}



"print.frailtyPenal" <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        if (x$type == "counting") {
            cat("\n      left truncated structure used")
        }
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" frailtyPenal failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- x$coef
    nvar <- length(x$coef)
 
    if (is.null(coef))
      {
        x$varH<-matrix(x$varH) 
        x$varHIH<-matrix(x$varHIH)
      }
    
    if (!is.null(coef)) 
      { 
        seH <- sqrt(diag(x$varH))[-1]
        seHIH <- sqrt(diag(x$varHIH))[-1]
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
        pchisq((coef/seH)^2, 1), digits - 1))
        cat("\n")
        if (!is.null(x$theta)) 
          {
          cat("  Shared Gamma Frailty model parameter estimates ", 
            "\n")
          cat("  using a Penalized Likelihood on the hazard function", 
            "\n")
          }
        else {
          cat("  Cox proportional hazards model parameter estimates ", 
            "\n")
          cat("  using a Penalized Likelihood on the hazard function", 
            "\n")
          }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "SE coef (H)", "SE coef (HIH)", "z", "p"))
        cat("\n")
        prmatrix(tmp)
        cat("\n")
      } 
        if (!is.null(x$theta)) {
          tetha <- x$theta
          temp <- diag(x$varH)[1]
          seH <- sqrt(((2 * (tetha^0.5))^2) * temp)
          temp <- diag(x$varHIH)[1]
          seHIH <- sqrt(((2 * (tetha^0.5))^2) * temp)
         cat("    Frailty parameter, Theta:", tetha, "(SE (H):", 
            seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        }
    cat(" \n")
    cat(paste("    penalized marginal log-likelihood =", round(x$logVerComPenal, 
        2)))
    cat("\n")
    cat("    n=", x$n)
    if (length(x$na.action)) 
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    else cat("\n")
    cat("    n events=", x$n.event, "\n")
    cat("    n groups=", x$groups, "\n")
    cat("    number of iterations: ", x$n.iter)
    cat("\n")
    cat("    Exact number of knots used: ", x$n.knots)
    cat("\n")
    invisible()
}



"summary.frailtyPenal"<-
 function(object,level=.95, len=6, d=2, lab="hr", ...)
    {
      x <- object
      if (!inherits(x, "frailtyPenal")) 
         stop("Invalid data")
      
      z<-abs(qnorm((1-level)/2))
      co <- x$coef
      se <- sqrt(diag(x$varH))[-1]
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
          }
      else
          lab <- dimnames(n)[[1]]
      
      mx <- max(nchar(lab)) + 1
      cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
      for(i in (1):dd[1]) {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }
      
}



"plot.frailtyPenal" <-
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


      } 
   }        

    return(invisible())
}


"lines.frailtyPenal"<-
function (x, type.plot="hazard", conf.bands=TRUE, ...) 
{

    plot.type <- charmatch(type.plot, c("hazard", "survival"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
       }

    if(plot.type==1)
       {
         if(conf.bands)
           {
             matlines(x$x1,x$lam,...) 
           }
         else
             lines(x$x1,x$lam[,1],...)   
       }   


    if(plot.type==2)
       {
         if(conf.bands)
           {
             matlines(x$x1,x$surv,...) 
           }
         else
             lines(x$x1,x$surv[,1],...)   
       }
  
  return(invisible())

}



