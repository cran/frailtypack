
"frailtyPenal" <-
function (formula, formula.terminalEvent, data, Frailty = FALSE, joint=FALSE, recurrentAG=FALSE, 
             cross.validation=FALSE, n.knots, kappa1 , kappa2, maxit=350)
 {

   
    if (missing(n.knots))
         stop("number of knots are required")       

    if (n.knots<4) n.knots<-4
    if (n.knots>20) n.knots<-20

    if (missing(kappa1))
         stop("smoothing parameter (kappa1) is required")       

    if (!missing(kappa2) & cross.validation)
        stop("The cross validation is not implemented for two strata")

    if (missing(kappa2) & joint)
        stop("smoothing parameter (kappa2) is required for the joint model")


    if (joint & cross.validation)
        stop("The cross validation is not implemented for the joint model")    


    call <- match.call()
    m <- match.call(expand = FALSE)
    m$formula.terminalEvent <- m$Frailty <- m$joint <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <-  m$... <- NULL
    special <- c("strata", "cluster", "subcluster", "terminal")
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
    subcluster <- attr(Terms, "specials")$subcluster
    terminalEvent <- attr(Terms, "specials")$terminal
   
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
 
    if(length(uni.cluster)>5000) 
     {
      stop("grouping variable must have less than 5000 groups
             \n please contact to the mantainer")   
     }
    

    if (length(subcluster)) {
       tempsub <- untangle.specials(Terms, "subcluster", 1:10)
        ordsub <- attr(Terms, "order")[tempsub$terms]
        if (any(ordsub > 1)) 
            stop("subcluster can not be used in an interaction")
        subcluster <- strata(m[, tempsub$vars], shortlabel = TRUE)
        dropx <- c(dropx,tempsub$terms)
        uni.subcluster<-unique(subcluster)

        if (joint)
         stop("joint model is not implemented for nested model")        


# fixed after Virginie's comment
        if (missing(kappa2))
         kappa2<-kappa1 

        if(length(uni.subcluster)==1) 
         {
           stop("subcluster variable must have more than 1 level")   
         }

        if(length(uni.subcluster)>5000) 
        {
         stop("sub-grouping variable must have less than 5000 groups
             \n please contact to the mantainer")   
        }

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
        if (!joint & !length(subcluster)) 
        kappa2<-0 
      }


    if (length(terminalEvent)) {
        tempterm <- untangle.specials(Terms, "terminal", 1:10)
        ord <- attr(Terms, "order")[tempterm$terms]
        if (any(ord > 1)) 
            stop("Terminal can not be used in an interaction")
        dropx <- c(dropx,tempterm$terms)
        terminal <- strata(m[, tempterm$vars], shortlabel = TRUE)
        terminal <- as.numeric(as.character(terminal))
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
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
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

 flush.console()
 ptm<-proc.time()
 cat("Be patient. The program is computing ... \n")


 if (!joint & !length(subcluster))  # shared model
  {
    ans <- .Fortran("frailpenal",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(cens.data),
                as.integer(uni.strat),
                as.integer(Frailty),
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
                as.integer(crossVal),
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
                as.integer(0),
                as.double(c(0,0)), 
                as.double(0),  PACKAGE = "frailtypack"
     )    
    
    if (noVar==1) nvar<-0

    np <- ans[[20]]
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans[[32]]
    fit$logVerComPenal <- ans[[24]]
    if (Frailty) {
        fit$theta <- (ans[[21]][np - nvar])^2
    }
    if (!Frailty) {
        fit$theta <- NULL
    }
    if (noVar==1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans[[21]][(np - nvar + 1):np]
       names(fit$coef) <- colnames(X)
     }
    
    temp1 <- matrix(ans[[22]], nrow = 50, ncol = 50)[1:np, 1:np]
    temp2 <- matrix(ans[[23]], nrow = 50, ncol = 50)[1:np, 1:np]
    fit$varH <- temp1[(np - nvar):np, (np - nvar):np]
    fit$varHIH <- temp2[(np - nvar):np, (np - nvar):np]
    fit$formula <- formula(Terms)
    fit$x1 <- ans[[25]]
    fit$lam <- matrix(ans[[26]], nrow = 99, ncol = 3)
    fit$surv <- matrix(ans[[27]], nrow = 99, ncol = 3)
    fit$x2 <- ans[[28]]
    fit$lam2 <- matrix(ans[[29]], nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans[[30]], nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots 
    fit$n.iter <- ans[[31]]
    fit$kappa <- ans[[34]]
    fit$DoF <- ans[[35]]
    fit$cross.Val<-cross.validation

    if(ans[[33]]==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    if(ans[[33]]==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")
    
    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-FALSE
    class(fit) <- "frailtyPenal"
 }

if (joint & !length(subcluster)) # joint model
 {
    if (!recurrentAG)
     {
      tt1.death<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
      tt0.death<-rep(0,length(tt1.death))
     }
    else
     {
     tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
     tt0.death<-rep(0,length(tt1.death))
     }
    

    Terms2 <- if (missing(data)) 
        terms(formula.terminalEvent, special)
    else terms(formula.terminalEvent, special, data = data)   
    ord2 <- attr(Terms2, "order")
    if (length(ord2) & any(ord2 != 1)) 
        stop("Interaction terms are not valid for terminal event formula")


    terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]

# Control sobre terminalEvent que debe ser 0-1 
    if (!all(terminalEvent%in%c(1,0))) 
        stop("terminal must contain a variable coded 0-1 and a non-factor variable")

    m2 <- match.call(expand = FALSE)

    m2$formula.terminalEvent <- m2$Frailty <- m2$joint <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <-  m2$... <- NULL

    m2$formula <- Terms2

    m2[[1]] <- as.name("model.frame")
    m2 <- eval(m2, sys.parent())


# Control sobre missing de recurrent para que cuadre el terminal
   match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]
   m2<-m2[match.noNA,]



##### New after Vigninie's comment and paper (for dealing with factors)

    newTerms2<-Terms2
    X2 <- model.matrix(newTerms2, m2)
    assign <- lapply(attrassign(X2, newTerms2)[-1], function(x) x - 1)
    if (ncol(X2) == 1) 
      {
         X2<-X2-1
      }
    else
      {
         X2 <- X2[, -1, drop = FALSE]
      } 

    nvar2<-ncol(X2) 

    if(nvar2>15)
       stop("maximum number of variables allowed for death are 15. 
             \n please contact to the mantainer")

    vardc.temp<-matrix(c(X2),nrow=nrow(X2),ncol=nvar2)
    
#####

   if(is.null(nrow(m2)))
        {
         if (length(m2) != nrow(m))
            stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")

        }
       else
        {

        if (nrow(m2) != nrow(m))
            stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")

        }


    if (!is.null(ncol(vardc.temp)))
     {
      vardc<-aggregate(vardc.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2] 
      for (i in 2:ncol(vardc.temp)) 
       {
        vardc.i<-aggregate(vardc.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2] 
        vardc<-cbind(vardc,vardc.i)
       }
     } 
    else
     {
      vardc<-aggregate(vardc.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]  
     }


# ojo, consideramos que puede que no sean las mismas variables para recurrent y para terminal event
# 1a version no, y tomaba nvar<-2*nvar
     
    nvarRec<-nvar
    if(is.null(nrow(vardc)))
     nvarEnd<-1
    else
     nvarEnd<-ncol(vardc)

    nvar<-nvarRec+nvarEnd

    ans <- .Fortran("frailpenalJoint",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(cens.data),
                as.integer(n.knots),
                k0=c(as.double(kappa1),as.double(kappa2)),
                as.double(tt0),
                as.double(tt1),
                as.integer(cens),
                as.integer(cluster),

                as.double(tt0.death),               
                as.double(tt1.death),
                as.integer(terminalEvent),

                as.integer(nvarRec),
                as.double(var),
                as.integer(nvarEnd),
                as.double(vardc),
                as.integer(AG),
                as.integer(noVar), 
                as.integer(maxit),
                as.integer(crossVal),

                np=as.integer(0),
                b=as.double(rep(0,50)),
                H=as.double(matrix(0,nrow=50,ncol=50)),
                HIH=as.double(matrix(0,nrow=50,ncol=50)),
                loglik=as.double(0),
                x1=as.double(rep(0,99)),
                lam=as.double(matrix(0,nrow=99,ncol=3)),
                surv=as.double(matrix(0,nrow=99,ncol=3)),
                x2=as.double(rep(0,99)),
                lam2=as.double(matrix(0,nrow=99,ncol=3)),
                surv2=as.double(matrix(0,nrow=99,ncol=3)),
                ni=as.integer(0),
                cpt=as.integer(0),
                cpt.dc=as.integer(0),
                ier=as.integer(0),
                PACKAGE = "frailtypack")    
    
    if (noVar==1) nvar<-0
    np <- ans$np
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans$cpt
    fit$n.deaths <- ans$cpt.dc
    fit$logVerComPenal <- ans$loglik



   
    fit$theta <- ans$b[np - nvar-1]^2
    fit$alpha <- ans$b[np - nvar]

    
    if (noVar==1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans$b[(np - nvar + 1):np]


       names(fit$coef) <- c(colnames(X), colnames(X2))
     }
    

    temp1 <- matrix(ans$H, nrow = 50, ncol = 50)[1:np, 1:np]
    temp2 <- matrix(ans$HIH, nrow = 50, ncol = 50)[1:np, 1:np]
    fit$nvar<-c(nvarRec,nvarEnd)
    fit$varH <- temp1[(np - (nvar) - 1):np, (np - (nvar) - 1):np]
    fit$varHIH <- temp2[(np - (nvar) - 1):np, (np - (nvar) - 1):np]
    fit$formula <- formula(Terms)
    fit$x1 <- ans$x1
    fit$lam <- matrix(ans$lam, nrow = 99, ncol = 3)
    fit$surv <- matrix(ans$surv, nrow = 99, ncol = 3)
    fit$x2 <- ans$x2
    fit$lam2 <- matrix(ans$lam2, nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans$surv2, nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.knots<-n.knots 
    fit$n.iter <- ans$ni
    fit$kappa <- ans$k0
    fit$cross.Val<-cross.validation

    if(ans$ier==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-FALSE
    class(fit) <- "jointPenal"

 }


if (length(subcluster))  # nested model
 {

    ans <- .Fortran("nested",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(length(uni.subcluster)),
                as.integer(cens.data),
                as.integer(uni.strat),
                as.integer(n.knots),
                as.double(kappa1),
                as.double(kappa2),
                as.double(tt0),
                as.double(tt1),
                as.integer(cens),
                as.integer(cluster),
                as.integer(subcluster),
                as.integer(nvar),
                as.double(strats),
                as.double(var),
                as.integer(AG),
                as.integer(noVar), 
                as.integer(maxit),
                as.integer(crossVal),  
                np=as.integer(0),
                b=as.double(rep(0,150)),
                H=as.double(matrix(0,nrow=50,ncol=50)),
                HIH=as.double(matrix(0,nrow=50,ncol=50)),
                loglik=as.double(0),
                x1=as.double(rep(0,99)),
                lam=as.double(matrix(0,nrow=99,ncol=3)),
                surv=as.double(matrix(0,nrow=99,ncol=3)),
                x2=as.double(rep(0,99)),
                lam2=as.double(matrix(0,nrow=99,ncol=3)),
                surv2=as.double(matrix(0,nrow=99,ncol=3)),
                ni=as.integer(0),
                cpt=as.integer(0),
                ier=as.integer(0),
                k0=as.double(c(0,0)), 
                ddl=as.double(0),  PACKAGE = "frailtypack")    


    if (noVar==1) nvar<-0

    np <- ans$np
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$subgroups <- length(uni.subcluster)
    fit$n.events <- ans$cpt
    fit$logVerComPenal <- ans$loglik
    
    fit$alpha<-ans$b[np-nvar-1]^2
    fit$eta<-ans$b[np-nvar]^2

    if (noVar==1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans$b[(np - nvar + 1):np]
       names(fit$coef) <- colnames(X)
     }
    

    temp1 <- matrix(ans$H, nrow = 50, ncol = 50)[1:np, 1:np]
    temp2 <- matrix(ans$HIH, nrow = 50, ncol = 50)[1:np, 1:np]
    fit$varH <- temp1[(np - nvar - 1):np, (np - nvar - 1):np]
    fit$varHIH <- temp2[(np - nvar - 1):np, (np - nvar - 1):np]
#    fit$varHIH <- fit$varH
    fit$formula <- formula(Terms)
    fit$x1 <- ans$x1
    fit$lam <- matrix(ans$lam, nrow = 99, ncol = 3)
    fit$surv <- matrix(ans$surv, nrow = 99, ncol = 3)
    fit$x2 <- ans$x2
    fit$lam2 <- matrix(ans$lam2, nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans$surv2, nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots 
    fit$n.iter <- ans$ni
    fit$kappa <- ans$k0
    fit$DoF <- ans$ddl
    fit$cross.Val<-cross.validation

    if(ans$ier==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")

    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-TRUE
    class(fit) <- "nestedPenal"
 }

 cost<-proc.time()-ptm
 cat("The program took", round(cost[3],2), "seconds \n")

 fit

    
}



