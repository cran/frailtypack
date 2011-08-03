
"frailtyPenal" <-
function (formula, formula.terminalEvent, data, Frailty = FALSE, joint=FALSE, recurrentAG=FALSE, 
             cross.validation=FALSE, n.knots, kappa1 , kappa2, maxit=350)
{

#AD:
	if (missing(formula))stop("The argument formula must be specified in any model")
	if(class(formula)!="formula")stop("The argument formula must be a formula")
#AD:   
	if (missing(n.knots))stop("number of knots are required")   
#AD:	 
	n.knots.temp <- n.knots	 
#AD
	if (n.knots<4) n.knots<-4
	if (n.knots>20) n.knots<-20

	if (missing(kappa1))stop("smoothing parameter (kappa1) is required")     
#AD:
	if (!missing(kappa2)) indic.Kappa2 <- 0
	if (missing(kappa2)) indic.Kappa2 <- 1
#AD:  

	if (!missing(kappa2) & cross.validation){
		stop("The cross validation is not implemented for two strata")
	}
	
	if (missing(kappa2) & joint){
		stop("smoothing parameter (kappa2) is required for the joint model")
	}


	if (joint & cross.validation){
		stop("The cross validation is not implemented for the joint model")  
	}  

	call <- match.call()
	
	m <- match.call(expand = FALSE)
    m$formula.terminalEvent <- m$Frailty <- m$joint <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <-  m$... <- NULL
  
    special <- c("strata", "cluster", "subcluster", "terminal")
		
	Terms <- if (missing(data)){ 
		terms(formula, special)
	}else{
		terms(formula, special, data = data)  
	} 
	
	ord <- attr(Terms, "order") ## longueur de ord=nbre de var.expli
		
	if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
#si pas vide tous si il ya au moins un qui vaut 1 on arrête
	
	m$formula <- Terms
	m[[1]] <- as.name("model.frame") ##m[[1]]=frailtypenal, il le remplace par model.frame en fait
		
#model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) + 
#as.factor(charlson) + sex + chemo + terminal(death), data = readmission)
		
	m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument
		
	if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
		
	Y <- model.extract(m, "response") # objet de type Surv =Time
		
	if (!inherits(Y, "Surv"))stop("Response must be a survival object") #test si c bien un objet de type "Surv"
		
	ll <- attr(Terms, "term.labels")#liste des variables explicatives
#cluster(id) as.factor(dukes) as.factor(charlson) sex chemo terminal(death) 
	
	mt <- attr(m, "terms") #m devient de class "formula" et "terms"
			
	X <- if (!is.empty.model(mt))model.matrix(mt, m, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables
			
	strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
	cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
	subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()
	terminalEvent <- attr(Terms, "specials")$terminal #nbre de var qui sont en fonction de terminal()
	
	dropx <- NULL

	if (length(cluster)){
		tempc <- untangle.specials(Terms, "cluster", 1:10)
		ord <- attr(Terms, "order")[tempc$terms]
		if (any(ord > 1))stop("Cluster can not be used in an interaction")
		cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
		dropx <- tempc$terms
		uni.cluster<-unique(cluster)
	}else{
		stop("grouping variable is needed")   
	}

	if(length(uni.cluster)==1){ 
		stop("grouping variable must have more than 1 level")   
	}

	if(length(uni.cluster)>1500 & missing(formula.terminalEvent)){ 
		stop("grouping variable must have less than 1500 groups
	\n please contact to the mantainer")   
	}

	if(length(uni.cluster)>15000 & !missing(formula.terminalEvent)) 
	{
		stop("grouping variable must have less than 15000 groups
	\n please contact to the mantainer")   
	}

		
	if (length(subcluster)){
		tempsub <- untangle.specials(Terms, "subcluster", 1:10)
		ordsub <- attr(Terms, "order")[tempsub$terms]		
		if (any(ordsub > 1))stop("subcluster can not be used in an interaction")
		subcluster <- strata(m[, tempsub$vars], shortlabel = TRUE)
		dropx <- c(dropx,tempsub$terms)
		uni.subcluster<-unique(subcluster)
		if (joint)stop("joint model is not implemented for nested model")        

		if (missing(kappa2))kappa2<-kappa1
	
		if(length(uni.subcluster)==1){        
			stop("subcluster variable must have more than 1 level")   
		}

		if(length(uni.subcluster)>5000){
			stop("sub-grouping variable must have less than 5000 groups
			\n please contact to the mantainer")   
		}

	}	
#AD:	
	if (length(cluster) == length(subcluster)){
		if (all(all.equal(cluster,subcluster)==T)){	
			stop("'Subgroup' variable and 'group' variable need to be different")
		}
	}	
#AD:	
	if (length(strats)){

		temp <- untangle.specials(Terms, "strata", 1)
		dropx <- c(dropx, temp$terms)
		if (length(temp$vars) == 1)strata.keep <- m[[temp$vars]]
		else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
		strats <- as.numeric(strata.keep)
		uni.strat<-length(unique(strats))

		if (missing(kappa1))stop("smoothing parameter (kappa1) is required") 

		if (uni.strat!=2)stop("maximum number of strata is 2")
	}else{
		uni.strat<-1
		strats <- rep(1,nrow(data))
		if (!joint & !length(subcluster)) 
		kappa2<-0 
	}
	
#AD: indicator of terminal()
	ind.terminal <- length(terminalEvent)
#AD:
	if (length(terminalEvent)){
	
		tempterm <- untangle.specials(Terms, "terminal", 1:10) 
#ici on comme terme tempterm$vars qui est le nom dans l'appel(ex;"terminal(death)"
#et tempterm$terms qui est la position de la variable dans l'appel, ici elle vient a la position 6
	
		ord <- attr(Terms, "order")[tempterm$terms] # ord[6]=1 ici dans notre exemple
		
		if (any(ord > 1))stop("Terminal can not be used in an interaction")
		dropx <- c(dropx,tempterm$terms) # vecteur de position
		terminal <- strata(m[, tempterm$vars], shortlabel = TRUE)
		terminal <- as.numeric(as.character(terminal))

	}
	
	type <- attr(Y, "type")

	if (type != "right" && type != "counting"){ 
		stop(paste("Cox model doesn't support \"", type, "\" survival data", 
		sep = ""))
	}

	if (type != "counting" && recurrentAG){
		stop("recurrentAG needs counting process formulation")
	}
	
	#drop contient les position liees au fonction() ic ex:cluster(id) et terminal(death)
	if (length(dropx)){ 
		newTerms <- Terms[-dropx]
	}else{
		newTerms <- Terms
	}
#newTerm vaut Terms - les variables dont les position sont dans drop
	X <- model.matrix(newTerms, m)  	
	
	assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
# assigne donne la position pour chaque variables
#ncol(X) : nombre de variable sans sans les fonction speciaux comme terminal()...+id
	

	if (ncol(X) == 1){ 
		X<-X-1
		noVar1 <- 1 
	}else{
		X <- X[, -1, drop = FALSE]
		noVar1 <- 0
	} 	
# on enleve ensuite la premiere colonne correspondant a id
	
	
	
	
	nvar<-ncol(X) #nvar==1 correspond a 2 situaions: 
# au cas ou on a aucune var explicative dans la partie rec, mais X=0
# cas ou on a 1seul var explicative, ici X est en general different de 0


	if(nvar>15){
		stop("maximum number of variables allowed are 15. 
			\n please contact to the mantainer")
	}
	
	
	var<-matrix(c(X),nrow=nrow(X),ncol=nvar) #matrix sans id et sans partie ex terminal(death)
	
	n<-nrow(X)    
	if(n>20000){     
		stop("number of observations must be less than 20000 
		\n please contact to the mantainer")   
	}
	


		#car Y est de class surv
	if (type=="right"){
		tt0 <- rep(0,n)
		tt1 <- Y[,1]
		cens <- Y[,2]
	}else{
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
	cat("\n")
	cat("Be patient. The program is computing ... \n")

	
	
		
#
# Begin SHARED MODEL
#

	


 if (!joint & !length(subcluster))  
  {
   if (sum(as.double(var))==0) nvar <- 0
   np <- as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + as.integer(Frailty)
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
                as.integer(noVar1), 
                as.integer(maxit),
                as.integer(crossVal),
                np=as.integer(np),
                as.double(rep(0,50)),
                as.double(matrix(0,nrow=np,ncol=np)),
                as.double(matrix(0,nrow=np,ncol=np)),
                as.double(0),
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
     
    if (noVar1 == 1) nvar<-0

    np <- ans[[20]]
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans[[33]]
    fit$logLikPenal <- ans[[24]]
    if (Frailty) {
        fit$theta <- (ans[[21]][np - nvar])^2
    }
    if (!Frailty) {
        fit$theta <- NULL
    }
    if (noVar1 == 1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans[[21]][(np - nvar + 1):np]
       names(fit$coef) <- colnames(X)
     }
#AD:modification des dimensions des tableaux
    temp1 <- matrix(ans[[22]], nrow = np, ncol = np)
    temp2 <- matrix(ans[[23]], nrow = np, ncol = np)
    fit$varH <- temp1[(np - nvar):np, (np - nvar):np]
    fit$varHIH <- temp2[(np - nvar):np, (np - nvar):np]
    fit$formula <- formula(Terms)
    fit$x1 <- ans[[26]]
    fit$lam <- matrix(ans[[27]], nrow = 99, ncol = 3)
    fit$surv <- matrix(ans[[28]], nrow = 99, ncol = 3)
    fit$x2 <- ans[[29]]
    fit$lam2 <- matrix(ans[[30]], nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans[[31]], nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots 
    fit$n.iter <- ans[[32]]
    fit$kappa <- ans[[35]]
    fit$DoF <- ans[[36]]
    fit$cross.Val<-cross.validation
#AD:
    fit$LCV <- ans[[25]]
    fit$npar <- np
    fit$nvar <- nvar
    fit$noVar1 <- noVar1
#AD:
    if(ans[[34]]==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    if(ans[[34]]==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")
 #AD:
	fit$n.knots.temp <- n.knots.temp
#AD   
    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-FALSE
    class(fit) <- "frailtyPenal"

 }  # End SHARED MODEL


#
# Begin JOINT MODEL
#
	
	if (joint & !length(subcluster)) 
	{

# Preparing data ...
#AD:
		if(Frailty =="FALSE"){
			stop("For joint frailty models, 'Frailty' must be equal to 'TRUE' ")
		}
#AD

		if (!recurrentAG)
		{
			tt1.death<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
			tt0.death<-rep(0,length(tt1.death))
     		}else{
			tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
			tt0.death<-rep(0,length(tt1.death))
		}

       
		Terms2 <- if (missing(data)){ 
		
			if (!missing(formula.terminalEvent))terms(formula.terminalEvent, special)
		}else{
			if (!missing(formula.terminalEvent))terms(formula.terminalEvent, special, data = data) 
		}
#AD:
		if (!missing(formula.terminalEvent)){
			ord2 <- attr(Terms2, "order")
		
			if (length(ord2) & any(ord2 != 1)){ 
				stop("Interaction terms are not valid for terminal event formula")
			}
		}
#AD:


#AD:Joint model needs "terminal()"
		if (ind.terminal){
			terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
		}else{
			stop(" Joint frailty model miss specified ")
		}
#AD:


# terminalEvent might be 0-1 
		if (!all(terminalEvent%in%c(1,0))){ 
			stop("terminal must contain a variable coded 0-1 and a non-factor variable")
		}
	
		m2 <- match.call(expand = FALSE)

## AD: modified 20 06 2011, for no covariates on terminal event part
		if (missing(formula.terminalEvent)){
			m2$Frailty <- m2$joint <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <-  m2$... <- NULL			
		}else{
			m2$formula.terminalEvent <- m2$Frailty <- m2$joint <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <-  m2$... <- NULL
		}

		
		m2$formula <- Terms2
		m2[[1]] <- as.name("model.frame")
		m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout
		

		match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA
		

		m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA
		
		if (!missing(formula.terminalEvent))newTerms2<-Terms2


		if (!missing(formula.terminalEvent)){
			X2 <- model.matrix(newTerms2, m2)
			assign <- lapply(attrassign(X2, newTerms2)[-1], function(x) x - 1)
			
			if (ncol(X2) == 1) 
			{
				X2<-X2-1
				noVar2 <- 1 
			}else{
				X2 <- X2[, -1, drop = FALSE]
				noVar2 <- 0 
			} 


			nvar2<-ncol(X2) 
		
			if(nvar2>15){
				stop("maximum number of variables allowed for death are 15. 
				\n please contact to the mantainer")
			}

			vardc.temp<-matrix(c(X2),nrow=nrow(X2),ncol=nvar2)


			if(is.null(nrow(m2)))
			{
				if (length(m2) != nrow(m)){
					stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
				}
			}else{
				
				if (nrow(m2) != nrow(m)){
					stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
				}
	
			}

			if (!is.null(ncol(vardc.temp))){
				vardc<-aggregate(vardc.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2] 
				
				if (ncol(vardc.temp)>1){
				
					for (i in 2:ncol(vardc.temp)){
						vardc.i<-aggregate(vardc.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2] 
						vardc<-cbind(vardc,vardc.i)
					} 
				}
			}else{
				vardc<-aggregate(vardc.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]
			}
		}else{
			noVar2 <- 1
			vardc<-0
		}

	nvarRec<-nvar

	if (!missing(formula.terminalEvent)){	
		if(is.null(nrow(vardc))){
			nvarEnd<-1
		}else{
			nvarEnd<-ncol(vardc)
		}
	}else{
		nvarEnd<-0
	}

	nvar<-nvarRec+nvarEnd 
	
	

# ... end preparing data 
#AD:
    effet <- 1
    indic_alpha <- 1
    nst=2
    np <- nst*(n.knots + 2) + nvarRec+nvarEnd +effet+indic_alpha

	if (all(all.equal(as.numeric(cens),terminal)==T)){	
		stop("'Recurrent event' variable and 'Terminal event' variable need to be different")    
	}

    
    ans <- .Fortran("joint",
                as.integer(n),
                as.integer(length(uni.cluster)), 

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
#AD:
                as.integer(noVar1), 
                as.integer(noVar2), 
#AD:
                as.integer(maxit),
                np=as.integer(np),		
                b=as.double(rep(0,np)),
                H=as.double(matrix(0,nrow=np,ncol=np)),
                HIH=as.double(matrix(0,nrow=np,ncol=np)),
                loglik=as.double(0),
		trace=as.double(0),
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
#AD:
  if(ans$ier==4 || ans$loglik==-1e9){
	cat("Problem in the loglikehood computation. The program stopped abnormally.\n")
	ans$H <- matrix(NA,nc=np,nr=np)
	ans$HIH <- matrix(NA,nc=np,nr=np)	
	ans$trace <- NA
	ans$x1 <- rep(NA,99)
	ans$lam <- matrix(NA,nr=99,nc=3)
	ans$surv <- matrix(NA,nr=99,nc=3)  
	ans$x2<- rep(NA,99)	
	ans$lam2 <- matrix(NA,nr=99,nc=3)  
	ans$surv2 <- matrix(NA,nr=99,nc=3)  
	ans$ni <- 0
	ans$cpt <- 0
	ans$cpt.dc <- 0
  }
#AD:	
#AD:    
    if (noVar1==1 & noVar2==1) nvar<-0
#AD:
    np <- ans$np
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans$cpt
    fit$n.deaths <- ans$cpt.dc
    fit$logLikPenal <- ans$loglik
#AD:
    fit$LCV <- ans$trace
#AD: 
    fit$theta <- ans$b[np - nvar-1]^2
    fit$alpha <- ans$b[np - nvar]
    fit$npar <- np

#AD:    
    if ((noVar1==1 & noVar2==1)) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans$b[(np - nvar + 1):np]

	if (missing(formula.terminalEvent)){
	   names(fit$coef) <- c(colnames(X)) 
	}else{
           names(fit$coef) <- c(colnames(X), colnames(X2))
        }
     }
    
#AD:
    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
    fit$nvar<-c(nvarRec,nvarEnd)
    fit$varH <- temp1[(np - (nvar) - 1):np, (np - (nvar) - 1):np]
    fit$varHIH <- temp2[(np - (nvar) - 1):np, (np - (nvar) - 1):np]
    fit$formula <- formula(Terms)
    fit$x1 <- ans$x1
    fit$lam <- matrix(ans$lam, nrow = 99, ncol = 3)
    fit$surv <- matrix(ans$surv, nrow = 99, ncol = 3)
    if (missing(formula.terminalEvent)){
    	fit$x2 <- NA
    }else{
        fit$x2 <- ans$x2
    }
    fit$lam2 <- matrix(ans$lam2, nrow = 99, ncol = 3)
    fit$surv2 <- matrix(ans$surv2, nrow = 99, ncol = 3)
    fit$type <- type
    fit$n.knots<-n.knots 
    fit$n.iter <- ans$ni
    fit$kappa <- ans$k0
    fit$cross.Val<-cross.validation
#AD:
    fit$noVar1 <- noVar1
    fit$noVar2 <- noVar2
    fit$nvarRec <- nvarRec
    fit$nvarEnd <- nvarEnd
#AD:

    if(ans$ier==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit){
        warning("model did not converge. Change the 'maxit' parameter") 
    }
#AD:
     fit$n.knots.temp <- n.knots.temp

#AD
    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-FALSE
    class(fit) <- "jointPenal"
 }  # End JOINT MODEL


#
# Begin NESTED MODEL
#

effet <- 1




if (length(subcluster))  
 {
    if (sum(as.double(var))==0) nvar <- 0
    np <- as.integer(uni.strat) * (as.integer(n.knots)+2) + as.integer(nvar) + 2*effet
    ans <- .Fortran("nested",
                as.integer(n),
                as.integer(length(uni.cluster)),
                as.integer(length(uni.subcluster)),
                as.integer(uni.strat),
                as.integer(n.knots),
                as.double(kappa1),
                as.double(kappa2),
                as.double(tt0),
                as.double(tt1),
                as.integer(cens),
                as.integer(cluster),
                as.integer(subcluster),
                nvar=as.integer(nvar),
                as.double(strats),
                as.double(var),
                as.integer(AG),
                as.integer(noVar1), 
                as.integer(maxit),
                as.integer(crossVal),  
                np=as.integer(np),
                b=as.double(rep(0,np)),
                H=as.double(matrix(0,nrow=np,ncol=np)),
                HIH=as.double(matrix(0,nrow=np,ncol=np)),
                loglik=as.double(0),
		trace=as.double(0),
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

#AD:
  if(ans$ier==4){
	cat("Problem in the loglikehood computation. The program stopped abnormally.\n")
	ans$H <- matrix(NA,nc=np,nr=np)
	ans$HIH <- matrix(NA,nc=np,nr=np)	
	ans$trace <- NA
	ans$x1 <- rep(NA,99)
	ans$lam <- matrix(NA,nr=99,nc=3)
	ans$surv <- matrix(NA,nr=99,nc=3)  
	ans$x2<- rep(NA,99)	
	ans$lam2 <- matrix(NA,nr=99,nc=3)  
	ans$surv2 <- matrix(NA,nr=99,nc=3)  
	ans$cpt <- 0
	ans$cpt.dc <- 0
  }
#AD:	
    nst <- as.integer(uni.strat)		
		
    if (noVar1 == 1) nvar<-0

    np <- ans$np
    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$subgroups <- length(uni.subcluster)
    fit$n.events <- ans$cpt
    fit$logLikPenal <- ans$loglik
    
    fit$alpha<-ans$b[np-nvar-1]^2
    fit$eta<-ans$b[np-nvar]^2

    if (noVar1 == 1) {
      fit$coef <- NULL
    } 
    else
     {
       fit$coef <- ans$b[(np - nvar + 1):np]
       names(fit$coef) <- colnames(X)
     }
    

    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
    fit$varH <- temp1[(np - nvar - 1):np, (np - nvar - 1):np]
    fit$varHIH <- temp2[(np - nvar - 1):np, (np - nvar - 1):np]
      
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
#AD:
    fit$nvar <- nvar
    fit$LCV <- ans$trace
    fit$npar <- np
    fit$nst <- nst
    fit$indic.Kappa2 <- indic.Kappa2 

#AD:

    if(ans$ier==-1)
        warning("matrix non-positive definite")

    if(fit$n.iter>maxit)
        warning("model did not converge. Change the 'maxit' parameter") 

    if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")
#AD:
     fit$n.knots.temp <- n.knots.temp

#AD
    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-TRUE
    class(fit) <- "nestedPenal"

 } # End NESTED MODEL

 cost<-proc.time()-ptm
 cat("The program took", round(cost[3],2), "seconds \n")
 fit

}



