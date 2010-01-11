
"additivePenal" <-
function (formula, data, correlation=FALSE, recurrentAG=FALSE, cross.validation=FALSE, n.knots, kappa1 , kappa2, maxit=350)
 {

    if (missing(n.knots))
         stop("number of knots are required")       

    if (n.knots<4) n.knots<-4
    if (n.knots>20) n.knots<-20

    if (missing(kappa1))
         stop("smoothing parameter (kappa1) is required")       

    if (!missing(kappa2) & cross.validation)
        stop("The cross validation is not implemented for two strata")

    call <- match.call()
    m <- match.call(expand = FALSE)
    m$correlation <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <-  m$... <- NULL
    special <- c("strata", "cluster", "slope")
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
    slope <- attr(Terms, "specials")$slope



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
 
 

    if (length(strats)) {
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
        uni.strat<-length(unique(strats))
      
        if (missing(kappa1))
         stop("smoothing parameter (kappa1) is required") 
    
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


    if (length(slope))
     {
        temps <- untangle.specials(Terms, "slope", 1:10)
        dropx <- c(dropx, temps$terms)
        ord <- attr(Terms, "order")[temps$terms]
        if (any(ord > 1)) 
            stop("'slope' can not be used in an interaction")
        aux<-temps$vars
        nnOK<-gsub(")","",gsub("slope\\(","",aux))
     }
    else
     {
      stop("interacting (slope) variable is needed")   
     }

    
    type <- attr(Y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("Cox model doesn't support \"", type, "\" survival data", 
            sep = ""))

    if (type != "counting" && recurrentAG)
       stop("recurrentAG needs counting process formulation")
    

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


    nn<-dimnames(X)[[2]]
    varInt<-c(1:length(nn))[nn==nnOK]


    nvar<-ncol(X) 

    if(nvar>15)
       stop("maximum number of variables allowed are 15. 
             \n please contact to the mantainer")

    var<-matrix(c(X),nrow=nrow(X),ncol=nvar)






    n<-nrow(X)    
    if(n>20000) 
     {
      stop("number of observations must be less than 20000 
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
    crossVal<-ifelse(cross.validation,0,1)

#    np<-uni.strat*n.knots+nvar+as.integer(correlation)+2*as.integer(frailty)+2
    np<-uni.strat*n.knots+nvar+as.integer(correlation)+2*1+2


    ptm<-proc.time()
    cat("Be patient. The program is computing ... \n")

    ans <- .Fortran("additive",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(cens.data),
                as.integer(uni.strat),
#                as.integer(frailty),
                as.integer(n.knots),
                as.double(kappa1),
                as.double(kappa2),
                as.double(tt0),
                as.double(tt1),
                as.integer(cens),
                as.integer(cluster),
                as.integer(nvar),
                as.integer(strats),
                as.double(var),
                as.integer(varInt), 
                as.integer(AG),
                as.integer(noVar), 
                as.integer(maxit),
                as.integer(crossVal),  
                as.integer(correlation),
                b=as.double(rep(0,150)),
                coef=as.double(rep(0,nvar)),   
                varcoef=as.double(rep(0,nvar)),
                varcoef2=as.double(rep(0,nvar)),
                rho=as.double(0),
                cov=as.double(0),
                varcov=as.double(0),
                varSigma2=as.double(c(0,0)),
                varTau2=as.double(c(0,0)),
                ni=as.integer(0),
                loglikpen=as.double(0),
                k0=as.double(c(0,0)),  
                x1=as.double(rep(0,99)),
                lam1=as.double(matrix(0,nrow=99,ncol=3)),
                surv1=as.double(matrix(0,nrow=99,ncol=3)),
                x2=as.double(rep(0,99)),
                lam2=as.double(matrix(0,nrow=99,ncol=3)),
                surv2=as.double(matrix(0,nrow=99,ncol=3)),
                ier=as.integer(0),
                ddl=as.double(0),
                PACKAGE = "frailtypack")   

    flush.console()
    cost<-proc.time()-ptm
    cat("The program took", round(cost[3],2), "seconds \n")

    if (noVar==1) nvar<-0

    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- sum(cens)  #coded 0: censure 1:event
    fit$logLikPenal <- ans$loglikpen
    fit$type <- type

    fit$b <- ans$b 

    if (noVar==1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans$coef
       names(fit$coef) <- colnames(X)
     }

    fit$varH <- ans$varcoef
    fit$varHIH <- ans$varcoef2

    fit$cov <- ans$cov
    fit$varcov <- ans$varcov

    fit$sigma2<-ans$b[np-nvar-1]^2
    fit$varSigma2<-ans$varSigma2
    fit$tau2<-ans$b[np-nvar]^2
    fit$varTau2<-ans$varTau2  
 
    fit$rho<-ans$rho
    fit$cov<-ans$cov
    fit$varcov<-ans$varcov    

    fit$formula <- formula(Terms)

    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots 
    fit$kappa <- ans$k0
    fit$n.iter <- ans$ni
    fit$cross.Val<-cross.validation
    fit$correlation<-correlation

    fit$x1 <- ans$x1
    fit$lam <- matrix(ans$lam1, nrow = 99, ncol = 3) 
    fit$surv <- matrix(ans$surv1, nrow = 99, ncol = 3)
    fit$x2 <- ans$x2
    fit$lam2 <- matrix(ans$lam2, nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans$surv2, nrow = 99, ncol = 3)


   
    fit$DoF <- ans$ddl

    if(ans$ier==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")


    class(fit) <- "additivePenal"
    fit
}


