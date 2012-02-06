
"frailtyPenal" <-
function (formula, formula.terminalEvent, data, Frailty = FALSE, joint=FALSE, recurrentAG=FALSE, 
             cross.validation=FALSE, n.knots, kappa1 , kappa2, maxit=350,hazard="Splines",nb.int1,nb.int2)
{
##### hazard specification ######
haztemp <- hazard
hazard <- strsplit(hazard,split="-")
hazard <- unlist(hazard)  
if(!(length(hazard) %in% c(1,2))){stop("Please check and revise the hazard argument according to the format specified in the help.")}

### longueur hazard = 1
if((all.equal(length(hazard),1)==T)==T){
   if(!(hazard %in% c("Weibull","Piecewise","Splines"))){
	stop("Only 'Weibull', 'Splines' or 'Piecewise' hazard can be specified in hazard argument.")
   }else{
	typeof <- switch(hazard,"Splines"=0,"Piecewise"=1,"Weibull"=2)	
	if(typeof %in% c(0,2)){
### Splines
		if (!(missing(nb.int1)) || !(missing(nb.int2))){
			stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int1' and 'nb.int2' arguments must be deleted.")
		}
		if (typeof == 0){
			size1 <- 100
			size2 <- 100
			equidistant <- 2	
			nbintervR <- 0
			nbintervDC <- 0
		}
### Weibull
		if (typeof == 2){
			equidistant <- 2
			nbintervR <- 0
			nbintervDC <- 0
			size1 <- 100
			size2 <- 100
		}
	}else{
		stop ("The hazard argument is incorrectly specified.Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")		
	}
   } 
}else{
### Picewise
      typeof <- 1
#### longueur hazard > 1
      if(!("Piecewise" %in% hazard)){
	stop("Only 'Piecewise' hazard can be specified in hazard argument in this case")
      }
      if(!(all(hazard %in% c("Piecewise","per","equi")))){
		stop ("The hazard argument is incorrectly specified.Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
      }else{
		if (!(haztemp %in% c("Piecewise-per","Piecewise-equi"))){
			stop ("The hazard argument is incorrectly specified. Please refer to the help file of frailtypack.")
		}
		equidistant <- switch(haztemp,"Piecewise-per"=0,"Piecewise-equi"=1)
      }
 }

#AD:
	if (missing(formula))stop("The argument formula must be specified in any model")
	if(class(formula)!="formula")stop("The argument formula must be a formula")
	
	if(typeof == 0){
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
	}else{
		if (!(missing(n.knots)) || !(missing(kappa1)) || !(missing(kappa2)) || !(missing(cross.validation))){
			stop("When parametric hazard function is specified, 'Kappa1', 'n.knots' and 'cross.validations' arguments must be deleted.")
		}
		n.knots <- 0
		kappa1 <- 0
		kappa2 <- 0
		crossVal <- 0
		
	}
	call <- match.call()
	
	m <- match.call(expand.dots = FALSE)
    m$formula.terminalEvent <- m$Frailty <- m$joint <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <-  m$... <- NULL
    

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


#=========================================================>


	mt <- attr(m, "terms") #m devient de class "formula" et "terms"
			
	X <- if (!is.empty.model(mt))model.matrix(mt, m, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables
			
#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
	strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
	cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
	subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

	if (length(subcluster)){
		ll <- ll[-grep("subcluster",ll)]
	}
	if (length(cluster)){
		ll <- ll[-grep("cluster",ll)]
	}
	if (length(strats)){
		ll <- ll[-grep("strata",ll)]
	}

	ind.place <- grep("factor",ll)

	vecteur <- NULL
	vecteur <- c(vecteur,ll[ind.place])

	mat.factor <- matrix(vecteur,ncol=1,nrow=length(vecteur))

 # Fonction servant a prendre les termes entre "as.factor"
	vec.factor <-apply(mat.factor,MARGIN=1,FUN=function(x){
	pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
	pos2 <- length(unlist(strsplit(x,split="")))-1
	return(substr(x,start=pos1,stop=pos2))})

	if(length(vec.factor) > 0){
		vect.fact <- attr(X,"dimnames")[[2]]

		vect.fact <- vect.fact[grep("factor",vect.fact)]
		vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
		pos2 <- length(unlist(strsplit(x,split="")))-2
		return(substr(x,start=pos1,stop=pos2))})		
		occur <- rep(0,length(vec.factor))
	

		for(i in 1:length(vec.factor)){
			occur[i] = sum(vec.factor[i] == vect.fact)
		}
	}

#=========================================================>	
	
	
	

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
	if(length(vec.factor) > 0){
#========================================>
		position <- unlist(assign,use.names=F)
	}
	
#========================================>

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
	
	var<-matrix(c(X),nrow=nrow(X),ncol=nvar) #matrix sans id et sans partie ex terminal(death)
	
	n<-nrow(X)    

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
	if (typeof == 0){
		crossVal<-ifelse(cross.validation,0,1)
	}	
	
	
	flush.console()
	ptm<-proc.time()
	cat("\n")
	cat("Be patient. The program is computing ... \n")
	
#=======================================>
#======= Construction du vecteur des indicatrice
	if(length(vec.factor) > 0){
#		ind.place <- ind.place -1
		k <- 0
		for(i in 1:length(vec.factor)){
			ind.place[i] <- ind.place[i]+k
				k <- k + occur[i]-1
		}
	}

#==================================
# Begin SHARED MODEL
#

 if (!joint & !length(subcluster))  
  {
   
        if(equidistant %in% c(0,1)){
		if (missing(nb.int1)) stop("Time interval 'nb.int1' is required")
		if (class(nb.int1) != "numeric") stop("The argument 'nb.int1' must be a numeric")	
		if ((nb.int1 < 1)) stop("Number of Time 'nb.int1' interval must be between 1 and 20")
		if (nb.int1 > 20){
			 nb.int1 <-20
			indic.nb.int1 <- 1 # equals 1 for nb.int1 > 20	 
		}else{
			indic.nb.int1 <- 0 # equals 1 for nb.int1 < 20
		}
		nbintervR <- nb.int1
		size1 <- 3*nbintervR
	}
	if ((typeof == 0) | (typeof == 2)) indic.nb.int1 <- 0 


	if (sum(as.double(var))==0) nvar <- 0

	np <- switch(as.character(typeof),
		"0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + as.integer(Frailty)),
		"1"=(as.integer(uni.strat) * nbintervR + nvar + as.integer(Frailty)),
		"2"=(as.integer(uni.strat) * 2 + nvar + as.integer(Frailty)))

		xSu1 <- rep(0,100)
		xSu2 <- rep(0,100)
		if (typeof==0){
			mt1 <- size1	
		}else{
			mt1 <- 100
		}
		size2 <- mt1

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
				b=as.double(rep(0,np)),
				as.double(matrix(0,nrow=np,ncol=np)),
				as.double(matrix(0,nrow=np,ncol=np)),
				as.double(0),
				LCV=as.double(rep(0,2)),
				as.double(rep(0,size1)),
				as.double(matrix(0,nrow=size1,ncol=3)),
				xSu1=as.double(xSu1),
				as.double(matrix(0,nrow=size2,ncol=3)),
				as.double(rep(0,size1)),
				as.double(matrix(0,nrow=size1,ncol=3)),
				xSu2=as.double(xSu2),
				as.double(matrix(0,nrow=size2,ncol=3)),
				as.integer(typeof),
				as.integer(equidistant),
				as.integer(nbintervR),
				as.integer(size1),
				as.integer(0),
				as.integer(0),
				as.integer(0),
				as.double(c(0,0)), 
				as.double(0),
				istop=as.integer(0),
				shape.weib=as.double(rep(0,2)),
				scale.weib=as.double(rep(0,2)),
				as.integer(mt1),
				zi=as.double(rep(0,(n.knots+6))),
				martingale.res=as.double(rep(0,as.integer(length(uni.cluster)))),
				martingaleCox=as.double(rep(0,n)),
				frailty.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.var=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.sd=as.double(rep(0,as.integer(length(uni.cluster)))),
				linear.pred=as.double(rep(0,n)),
				time=as.double(rep(0,(nbintervR+1))),
				PACKAGE = "frailtypack")    
#AD:

    if (ans$istop == 4){
	 warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")    
     }
    
    if (ans$istop == 2){
         warning("Model did not converge. Change the 'maxit' parameter")
    }
    if (ans$istop == 3){
         warning("Matrix non-positive definite.")
    }
    
#AD:  

    if (noVar1 == 1) nvar<-0
    
    np <- ans[[20]]
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans[[39]]
    if(as.character(typeof)=="0"){    
        fit$logLikPenal <- ans[[24]]
    }else{
        fit$logLik <- ans[[24]]
    }
    
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

    fit$varH <- temp1[(np - nvar + 1):np, (np - nvar + 1):np]
    fit$varHIH <- temp2[(np - nvar + 1):np, (np - nvar + 1):np]
    if (Frailty) fit$varTheta <- c(temp1[(np - nvar),(np - nvar)],temp2[(np - nvar),(np - nvar)])
   

    fit$formula <- formula(Terms)

    fit$x1 <- ans[[26]]
    fit$lam <- matrix(ans[[27]], nrow = size1, ncol = 3)
    fit$x2 <- ans[[30]]
    fit$lam2 <- matrix(ans[[31]], nrow = size1, ncol = 3)
    
    fit$surv <- matrix(ans[[29]], nrow = size2, ncol = 3)
    fit$surv2 <- matrix(ans[[33]], nrow = size2, ncol = 3)

    fit$xSu1 <- ans$xSu1
    fit$xSu2 <- ans$xSu2

    fit$type <- type
    fit$n.strat <- uni.strat
    
     
    fit$n.iter <- ans[[38]]

    if (typeof == 0){
    	fit$n.knots<-n.knots
	fit$kappa <- ans[[41]]    
	fit$DoF <- ans[[42]]
	fit$cross.Val<-cross.validation
	fit$n.knots.temp <- n.knots.temp
	fit$zi <- ans$zi
    }
	if(typeof == 1) fit$time <- ans$time
#AD:
    
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
    fit$npar <- np
    fit$nvar <- nvar
    fit$noVar1 <- noVar1
    fit$indic.nb.int1 <- indic.nb.int1
#AD:
 
    if(ans[[40]]==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")

    fit$typeof <- typeof
    fit$equidistant <- equidistant
    fit$nbintervR <- nbintervR
    fit$istop <- ans$istop
    
    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib
    
    if(Frailty){  
	fit$martingale.res <- ans$martingale.res
	fit$frailty.pred <- ans$frailty.pred
	fit$frailty.var <- ans$frailty.var	
	fit$frailty.sd <- ans$frailty.sd
	
    }else{

	fit$martingaleCox <- ans$martingaleCox
    }
	fit$linear.pred <- ans$linear.pred 
#
#========================= Test de Wald pour shared

	
	if(length(vec.factor) > 0){
		Beta <- ans[[21]][(np - nvar + 1):np]
		VarBeta <- fit$varH#[2:(nvar+1),2:(nvar+1)] 
		nfactor <- length(vec.factor)
		p.wald <- rep(0,nfactor)
		
		fit$global_chisq <- waldtest(N=nvar,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta)
		fit$dof_chisq <- occur
		fit$global_chisq.test <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factor)){
			p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
		}
		fit$p.global_chisq <- p.wald
		fit$names.factor <- vec.factor 

		
	}else{
		fit$global_chisq.test <- 0
	}
	
#===============================================	


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
	
		m2 <- match.call(expand.dots = FALSE)

## AD: modified 20 06 2011, for no covariates on terminal event part
		if (missing(formula.terminalEvent)){
			m2$Frailty <- m2$joint <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <- m2$nb.int2 <-  m2$... <- NULL			
		}else{
			m2$formula.terminalEvent <- m2$Frailty <- m2$joint <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <- m2$nb.int2 <-  m2$... <- NULL
		}

		
		m2$formula <- Terms2
		m2[[1]] <- as.name("model.frame")
		m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout
		

		match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA
		

		m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA
		
		if (!missing(formula.terminalEvent))newTerms2<-Terms2

#=========================================================>
                lldc <- attr(newTerms2,"term.labels")
		ind.placedc <- grep("factor",lldc)
		vecteur <- NULL
		vecteur <- c(vecteur,lldc[ind.placedc])
		mat.factor <- matrix(vecteur,ncol=1,nrow=length(vecteur))

 # Fonction servant a prendre les termes entre "as.factor"
		vec.factordc <-apply(mat.factor,MARGIN=1,FUN=function(x){
		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
		pos2 <- length(unlist(strsplit(x,split="")))-1
		return(substr(x,start=pos1,stop=pos2))})	
	       
#=========================================================>
		if (!missing(formula.terminalEvent)){
			X2 <- model.matrix(newTerms2, m2)
#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
			if(length(vec.factordc) > 0){
				vect.fact <- attr(X2,"dimnames")[[2]]
				occurdc <- rep(0,length(vec.factordc))
	
				for(i in 1:length(vec.factordc)){
					occurdc[i] = length(grep(vec.factordc[i],vect.fact))
				}
			}
#=========================================================>
			assign <- lapply(attrassign(X2, newTerms2)[-1], function(x) x - 1)
#========================================>
			if(length(vec.factordc) > 0){
				positiondc <- unlist(assign,use.names=F)
			}
#========================================>
			if (ncol(X2) == 1) 
			{
				X2<-X2-1
				noVar2 <- 1 
			}else{
				X2 <- X2[, -1, drop = FALSE]
				noVar2 <- 0 
			} 


			nvar2<-ncol(X2) 
		
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
#=======================================>
#======= Construction du vecteur des indicatrice

		if(length(vec.factordc) > 0){
			k <- 0
			for(i in 1:length(vec.factordc)){
				ind.placedc[i] <- ind.placedc[i]+k
					k <- k + occurdc[i]-1
			}
		}

#==================================
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
	
        if(equidistant %in% c(0,1)){
		if (missing(nb.int1)) stop("Recurrent Time interval 'nb.int1' is required")
		if (class(nb.int1) != "numeric") stop("The argument 'nb.int1' must be a numeric")
		if (nb.int1 < 1) stop("Number of Time interval 'nb.int1' must be between 1 and 20")
		if (missing(nb.int2)) stop("Death Time intervale 'nb.int2' is required")
		if (class(nb.int2) != "numeric") stop("The argument 'nb.int2' must be a numeric")
		if (nb.int2 < 1) stop("Number of Time interval 'nb.int2' must be between 1 and 20")
		
		if (nb.int1 > 20){
			 nb.int1 <-20
			indic.nb.int1 <- 1 # equals 1 for nb.int1 > 20	 
		}else{
			indic.nb.int1 <- 0 # equals 1 for nb.int1 < 20
		}
		
		if (nb.int2 > 20){
			 nb.int2 <-20
			indic.nb.int2 <- 1 # equals 1 for nb.int1 > 20	 
		}else{
			indic.nb.int2 <- 0 # equals 1 for nb.int1 < 20
		}
		
		nbintervR <- nb.int1
		size1 <- 3*nbintervR	
		nbintervDC <- nb.int2
		size2 <- 3*nbintervDC		
	}
	if ((typeof == 0) | (typeof == 2)){
		indic.nb.int1 <- 0
		indic.nb.int2 <- 0
	}	
	np <- switch(as.character(typeof),
		"0"=(nst * (n.knots + 2) + nvarRec + nvarEnd + effet + indic_alpha),

		"1"=(nbintervR + nbintervDC + nvarRec + nvarEnd + effet + indic_alpha),

		"2"=(2*nst + nvar + effet + indic_alpha))


	if (all(all.equal(as.numeric(cens),terminal)==T)){	
		stop("'Recurrent event' variable and 'Terminal event' variable need to be different")    
	}

	xSu1 <- rep(0,100)
	xSu2 <- rep(0,100)
	if (typeof==0){
		mt11 <- size1	
		mt12 <- size2
	}else{
		mt11 <- 100
		mt12 <- 100
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
		LCV=as.double(rep(0,2)),
                x1=as.double(rep(0,size1)),
                lam=as.double(matrix(0,nrow=size1,ncol=3)),
                xSu1=as.double(xSu1),
                surv=as.double(matrix(0,nrow=mt11,ncol=3)),
                x2=as.double(rep(0,size2)),
                lam2=as.double(matrix(0,nrow=size2,ncol=3)),
		xSu2=as.double(xSu2),
                surv2=as.double(matrix(0,nrow=mt12,ncol=3)),
#
		as.integer(typeof),
		as.integer(equidistant),
		as.integer(nbintervR),
		as.integer(nbintervDC),
		as.integer(size1),
		as.integer(size2),
#
                ni=as.integer(0),
                cpt=as.integer(0),
                cpt.dc=as.integer(0),
                ier=as.integer(0),
		istop=as.integer(0),
		shape.weib=as.double(rep(0,2)),
		scale.weib=as.double(rep(0,2)),
		as.integer(mt11),
		as.integer(mt12),
		
		martingale.res=as.double(rep(0,as.integer(length(uni.cluster)))),
		martingaledc.res=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.var=as.double(rep(0,as.integer(length(uni.cluster)))),
		linear.pred=as.double(rep(0,n)),
		lineardc.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
		zi=as.double(rep(0,(n.knots+6))),
		time=as.double(rep(0,(nbintervR+1))),
		timedc=as.double(rep(0,(nbintervDC+1))),
                PACKAGE = "frailtypack")  	
		


    if (ans$istop == 4){
	 warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")    
     }
    
    if (ans$istop == 2){
         warning("Model did not converge. Change the 'maxit' parameter")
    }
    if (ans$istop == 3){
         warning("Matrix non-positive definite.")
    }

#AD:    
    if (noVar1==1 & noVar2==1) nvar<-0
#AD:
    
    np <- ans$np
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- ans$cpt
    fit$n.deaths <- ans$cpt.dc

    if(as.character(typeof)=="0"){    
        fit$logLikPenal <- ans$loglik
    }else{
        fit$logLik <- ans$loglik
    }    
#AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
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
    fit$lam <- matrix(ans$lam, nrow = size1, ncol = 3)
    fit$surv <- matrix(ans$surv, nrow = mt11, ncol = 3)
    if (missing(formula.terminalEvent)){
    	fit$x2 <- NA
    }else{
        fit$x2 <- ans$x2
    }
    fit$lam2 <- matrix(ans$lam2, nrow = size2, ncol = 3)
    fit$surv2 <- matrix(ans$surv2, nrow = mt12, ncol = 3)

    fit$xSu1 <- ans$xSu1
    fit$xSu2 <- ans$xSu2

    fit$type <- type
    fit$n.iter <- ans$ni
    fit$typeof <- typeof
    if (typeof == 0){  
	fit$n.knots<-n.knots      
	fit$kappa <- ans$k0
	fit$cross.Val<-cross.validation
	fit$n.knots.temp <- n.knots.temp	
	fit$zi <- ans$zi
    }	
    if(typeof == 1){
	fit$time <- ans$time
	fit$timedc <- ans$timedc
    }
#AD:
    fit$noVar1 <- noVar1
    fit$noVar2 <- noVar2
    fit$nbintervR <- nbintervR
    fit$nbintervDC <- nbintervDC
    fit$nvarRec <- nvarRec
    fit$nvarEnd <- nvarEnd
    fit$istop <- ans$istop
    fit$indic.nb.int1 <- indic.nb.int1
    fit$indic.nb.int2 <- indic.nb.int2
    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib
#AD:

    if (Frailty){    
		
	fit$martingale.res <- ans$martingale.res
	fit$martingaledeath.res <- ans$martingaledc.res
	
	fit$frailty.pred <- ans$frailty.pred
	fit$frailty.var <- ans$frailty.var
	
	fit$linear.pred <- ans$linear.pred  
	fit$lineardeath.pred <- ans$lineardc.pred
    }    
 


#================================> For the reccurrent
#========================= Test de Wald pour shared
	

	if(length(vec.factor) > 0){
		Beta <- ans$b[(np - nvar + 1):np]
		VarBeta <- diag(diag(fit$varH)[-c(1,2)])
		nfactor <- length(vec.factor)
		p.wald <- rep(0,nfactor)
		ntot <- nvarEnd + nvarRec
		fit$global_chisq <- waldtest(N=nvarRec,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarEnd,Ntot=ntot)
		fit$dof_chisq <- occur
		fit$global_chisq.test <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factor)){
			p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
		}
		fit$p.global_chisq <- p.wald
		fit$names.factor <- vec.factor 

		
	}else{
		fit$global_chisq.test <- 0
	}

#================================> For the death
#========================= Test de Wald pour shared

	if(length(vec.factordc) > 0){
		Beta <- ans$b[(np - nvar + 1):np]
		VarBeta <- diag(diag(fit$varH)[-c(1,2)]) 
		nfactor <- length(vec.factordc)
		p.walddc <- rep(0,nfactor)
		ntot <- nvarEnd + nvarRec
		fit$global_chisq_d <- waldtest(N=nvarEnd,nfact=nfactor,place=ind.placedc,modality=occurdc,b=Beta,Varb=VarBeta,Lfirts=nvarRec,Ntot=ntot)
		fit$dof_chisq_d <- occurdc
		fit$global_chisq.test_d <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factordc)){
			p.walddc[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurdc[i]), 3)
		}
		fit$p.global_chisq_d <- p.walddc
		fit$names.factordc <- vec.factordc 

		
	}else{
		fit$global_chisq.test_d <- 0
	}	
  
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
         if(equidistant %in% c(0,1)){
		if (missing(nb.int1)) stop("Time interval 'nb.int1' is required")
		if (class(nb.int1) != "numeric") stop("The argument 'nb.int1' must be a numeric")	
		if (nb.int1 < 1) stop("Number of Time 'nb.int2' interval must be between 1 and 20")
		if (nb.int1 > 20){
			 nb.int1 <-20
			indic.nb.int1 <- 1 # equals 1 for nb.int1 > 20	 
		}else{
			indic.nb.int1 <- 0 # equals 1 for nb.int1 < 20
		}
		nbintervR <- nb.int1
		size1 <- 3*nbintervR
	}
	 if ((typeof == 0) | (typeof == 2)) indic.nb.int1 <- 0
	if (sum(as.double(var))==0) nvar <- 0
   
	np <- switch(as.character(typeof),
		"0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + 2 * as.integer(Frailty)),

		"1"=(as.integer(uni.strat) * nbintervR + nvar + 2 * as.integer(Frailty)),

		"2"=(as.integer(uni.strat) * 2 + nvar + 2 * as.integer(Frailty)))   
    
		xSu1 <- rep(0,100)
		xSu2 <- rep(0,100)
		if (typeof==0){
			mt1 <- size1	
		}else{
			mt1 <- 100
		}
		size2 <- mt1
	



########### group and subgroup
	grpe <- function(g){

		grp <- unique(g)

		res <- rep(0,length(grp))	

		for(i in 1:length(res)){
			res[i] = sum(grp[i]==g)
		}
		return(res)
	}

	grp <- grpe(as.integer(cluster))

	subgrpe <- function(g,sg){

		j <- 0
		k <- 0
		res <- rep(0,length(g))

		for(i in 1:length(g)){
			k <- k + g[i]
			j <- j + 1
			temp <- sg[j:k]
			res[i] <- length(grpe(temp))
			j <- k
		}
		return(res)
	}

	subgbyg <- subgrpe(grp,as.integer(subcluster))
	maxng <- max(subgbyg)
	ngg <- length(uni.cluster)
#	cat("nombre de sujet par groupe\n")
#	print(grp)
#	cat("nombre de sous-groupe par groupe\n")
#	print(subgbyg)	
#### group and subgroup

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
                as.integer(nvar),
                as.double(strats),
                as.double(var),
                as.integer(AG),
                as.integer(noVar1), 
                as.integer(maxit),
                as.integer(crossVal),  
                as.integer(np),
                as.integer(maxng),
                b=as.double(rep(0,np)),
                H=as.double(matrix(0,nrow=np,ncol=np)),
                HIH=as.double(matrix(0,nrow=np,ncol=np)),
                loglik=as.double(0),
		LCV=as.double(rep(0,2)),
                x1=as.double(rep(0,size1)),
                lam=as.double(matrix(0,nrow=size1,ncol=3)),
		xSu1=as.double(xSu1),
                surv=as.double(matrix(0,nrow=size2,ncol=3)),
                x2=as.double(rep(0,size1)),
                lam2=as.double(matrix(0,nrow=size1,ncol=3)),
		xSu2=as.double(xSu2),
                surv2=as.double(matrix(0,nrow=size2,ncol=3)),
		
		as.integer(typeof),
		as.integer(equidistant),
		as.integer(nbintervR),
		as.integer(size1),
			
                ni=as.integer(0),
                cpt=as.integer(0),
                ier=as.integer(0),
                k0=as.double(c(0,0)), 
                ddl=as.double(0),
		istop=as.integer(0),
		shape.weib=as.double(rep(0,2)),
		scale.weib=as.double(rep(0,2)),
		as.integer(mt1),
		zi=as.double(rep(0,(n.knots+6))),
		time=as.double(rep(0,(nbintervR+1))),

		martingale.res=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.pred.group=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.pred.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
		frailty.var.group=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.var.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
		frailty.sd.group=as.double(rep(0,as.integer(length(uni.cluster)))),
		frailty.sd.subgroup=as.double(matrix(0,nrow=ngg,ncol=maxng)),
		linear.pred=as.double(rep(0,n)),
		PACKAGE = "frailtypack")    

    if (ans$istop == 4){
	 warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")    
     }
    
    if (ans$istop == 2){
         warning("Model did not converge. Change the 'maxit' parameter")
    }
    if (ans$istop == 3){
         warning("Matrix non-positive definite.")
    }

    nst <- as.integer(uni.strat)		
		
    if (noVar1 == 1) nvar<-0
   
    np <- np
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$subgroups <- length(uni.subcluster)
    fit$n.events <- ans$cpt 
    if(as.character(typeof)=="0"){    
        fit$logLikPenal <- ans$loglik
    }else{
        fit$logLik <- ans$loglik
    }     
     
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
    fit$x2 <- ans$x2
    fit$lam <- matrix(ans$lam, nrow = size1, ncol = 3)
    fit$lam2 <- matrix(ans$lam2, nrow = size1, ncol = 3)

    fit$surv <- matrix(ans$surv, nrow = size2, ncol = 3) 
    fit$surv2 <- matrix(ans$surv2, nrow = size2, ncol = 3)

    fit$xSu1 <- ans$xSu1
    fit$xSu2 <- ans$xSu2

    fit$type <- type
    fit$n.strat <- uni.strat
    fit$n.iter <- ans$ni
    fit$typeof <- typeof
    fit$noVar1 <- noVar1
    
    if (typeof == 0){
    	fit$n.knots<-n.knots
    	fit$kappa <- ans$k0
    	fit$DoF <- ans$ddl
    	fit$cross.Val<-cross.validation
	fit$zi <- ans$zi
    }
    if(typeof == 1)fit$time <- ans$time
#AD:
    fit$nbintervR <- nbintervR
    fit$nvar <- nvar
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]    
    fit$npar <- np
    fit$nst <- nst
    if (typeof == 0){
    	fit$indic.Kappa2 <- indic.Kappa2 
	fit$n.knots.temp <- n.knots.temp
    }
    fit$indic.nb.int1 <- indic.nb.int1
    fit$istop <- ans$istop
    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib
 
 #   if (Frailty){    
	fit$martingale.res <- ans$martingale.res
	fit$frailty.pred.group <- ans$frailty.pred.group

	nom1 <- paste("g_",c(1:ngg),sep="")
	nom2 <- paste("sub_g",c(1:maxng))

	frailty.pred.subgroup <- as.data.frame(matrix(round(ans$frailty.pred.subgroup,6),ncol=maxng))
	rownames(frailty.pred.subgroup) <- nom1
	colnames(frailty.pred.subgroup) <- nom2
	if(sum(which(subgbyg < max(subgbyg)))>0)frailty.pred.subgroup[which(subgbyg < max(subgbyg)),(subgbyg[which(subgbyg < max(subgbyg))]+1):max(subgbyg)] <- "."
	fit$frailty.pred.subgroup <- frailty.pred.subgroup

#	fit$frailty.var.group <- ans$frailty.var.group
	
# 	frailty.var.subgroup <- as.data.frame(matrix(round(ans$frailty.var.subgroup,6),nc=maxng))
# 	rownames(frailty.var.subgroup) <- nom1
# 	colnames(frailty.var.subgroup) <- nom2
# 
# 	if(sum(which(subgbyg < max(subgbyg)))>0)frailty.var.subgroup[which(subgbyg < max(subgbyg)),(subgbyg[which(subgbyg < max(subgbyg))]+1):max(subgbyg)] <- "."
	
#	fit$frailty.var.subgroup <- frailty.var.subgroup
	
#	fit$frailty.sd.group <- ans$frailty.sd.group
	
#	frailty.sd.subgroup <- as.data.frame(matrix(round(ans$frailty.sd.subgroup,6),nc=maxng))
#	rownames(frailty.sd.subgroup) <- nom1
#	colnames(frailty.sd.subgroup) <- nom2
#	if(sum(which(subgbyg < max(subgbyg)))>0)frailty.sd.subgroup[which(subgbyg < max(subgbyg)),(subgbyg[which(subgbyg < max(subgbyg))]+1):max(subgbyg)] <- "."
	
#	fit$frailty.sd.subgroup <- frailty.sd.subgroup


	fit$linear.pred <- ans$linear.pred  
	fit$subgbyg <- subgbyg
 #   }      
    if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")




#================================> For the reccurrent

#========================= Test de Wald pour shared



	if(length(vec.factor) > 0){
		Beta <- ans$b[(np - nvar + 1):np]
		VarBeta <- fit$varH[2:(nvar+1),2:(nvar+1)] 

		nfactor <- length(vec.factor)
		p.wald <- rep(0,nfactor)
		
		fit$global_chisq <- waldtest(N=nvar,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta)
		fit$dof_chisq <- occur
		fit$global_chisq.test <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factor)){
			p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
		}
		fit$p.global_chisq <- p.wald
		fit$names.factor <- vec.factor 

		
	}else{
		fit$global_chisq.test <- 0
	}


    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-TRUE
    class(fit) <- "nestedPenal"

 } # End NESTED MODEL

 cost<-proc.time()-ptm
 cat("The program took", round(cost[3],2), "seconds \n")
 fit

}



