
"frailtyPenal" <-
function (formula, formula.terminalEvent, data, recurrentAG=FALSE, cross.validation=FALSE, n.knots, kappa1 , kappa2,
             maxit=350, hazard="Splines", nb.int1, nb.int2, RandDist="Gamma", betaknots=1, betaorder=3, 
             init.B, init.Theta, init.Alpha, Alpha, LIMparam=1e-4, LIMlogl=1e-4, LIMderiv=1e-4, print.times=TRUE)
{

# al suppression de l'argument joint
if (!missing(formula.terminalEvent)) joint <- TRUE
else joint <- FALSE
if ((!missing(Alpha) | !missing(init.Alpha)) & !joint) stop("init.Alpha and Alpha parameters belong to joint frailty model")

#ad 15/02/12 :add Audrey
m2 <- match.call()
m2$formula <- m2$formula.terminalEvent <- m2$recurrentAG <- m2$cross.validation <- m2$n.knots <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <-m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$... <- NULL
Names.data <- m2$data

#### Betaknots et betaorder ####
if (betaknots > 10) stop("Number of knots for beta(t) greater than 10 is useless, please choose a number between 0 and 10 (3 is optimal)")
if ((betaorder == 0)|(betaorder > 4)) stop("B-splines order for beta(t) must be a number between 1 and 4 (3 is optimal)")

#### Frailty distribution specification ####
if (!(RandDist %in% c("Gamma","LogN"))) { stop("Only 'Gamma' and 'LogN' distributions for frailties are allowed") }
logNormal <- switch(RandDist,"Gamma"=0,"LogN"=1)

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
		stop ("The hazard argument is incorrectly specified. Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
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
			stop("When parametric hazard function is specified, 'Kappa1', 'n.knots' and 'cross.validation' arguments must be deleted.")
		}
		n.knots <- 0
		kappa1 <- 0
		kappa2 <- 0
		crossVal <- 0
		
	}
	call <- match.call()
	
	m <- match.call(expand.dots = FALSE) # recupere l'instruction de l'utilisateur

	m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <-  m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL


	special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")
	
	Terms <- if (missing(data)){ 
		terms(formula, special)
	}else{
		terms(formula, special, data = data)
	}
	
	ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli
	
	if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
#si pas vide tous si il ya au moins un qui vaut 1 on arrete
	
	m$formula <- Terms
	
	
	m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait

#model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) +
#as.factor(charlson) + sex + chemo + terminal(death), data = readmission)
	
	m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument

	cluster <- attr(Terms, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()

	# al 13/02/14 : suppression de l'argument Frailty
	if (length(cluster)) Frailty <- TRUE
	else Frailty <- FALSE

	subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

	# booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
	classofY <- attr(model.extract(m, "response"),"class")
	# attention le package pec rajoute un element dans l'attribut "class" des objets de survie
	if (length(classofY)>1) classofY <- classofY[2]

	typeofY <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..
	
#Al : tri du jeu de donnees par cluster croissant
	if (length(cluster)){
		tempc <- untangle.specials(Terms, "cluster", 1:10)
		ord <- attr(Terms, "order")[tempc$terms]
		if (any(ord > 1))stop("Cluster can not be used in an interaction")
		m <- m[order(m[,tempc$vars]),] # soit que des nombres, soit des caracteres
		ordre <- as.integer(row.names(m)) # recupere l'ordre du data set
		cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
		uni.cluster <- unique(cluster)
	}

# verification de la sutructure nested si besoin
	if (length(subcluster) && Frailty == TRUE){
		tempsub <- untangle.specials(Terms, "subcluster", 1:10)
		ordsub <- attr(Terms, "order")[tempsub$terms]
		if (any(ordsub > 1))stop("subcluster can not be used in an interaction")

		if (any(ifelse(apply(ifelse(table(m[,tempsub$vars],m[,tempc$vars])>0,1,0),1,sum)==1,FALSE,TRUE))){
			stop("nested structure is necessary to fit a nested model")
		}

		# tri par ordre croissant de subcluster a l'interieur des clusters
		m <- m[order(m[,tempc$vars],m[,tempsub$vars]),]
		subcluster <- strata(m[, tempsub$vars], shortlabel = TRUE)
		ordre <- as.integer(row.names(m))

		subcluster <- as.integer(subcluster) # a determiner si il y en a besoin
		curr <- subcluster[1]
		subcluster[1] <- 1
		for (i in 2:length(subcluster)) {
			if (subcluster[i] == curr) { subcluster[i] <- subcluster[i-1] }
			else {
				curr <- subcluster[i]
				subcluster[i] <- subcluster[i-1] + 1
			}
		}
	}
#Al

	if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0
		
	Y <- model.extract(m, "response") # objet de type Surv =Time

	if (classofY == "SurvIC") intcens <- TRUE # booleen censure par intervalle
	else intcens <- FALSE

	if (intcens == TRUE) {
		if (classofY != "SurvIC") stop("When interval censoring, must use the SurvIC fonction")
	} else {
		#if (!inherits(Y, "Surv")) stop("Response must be a survival object") #test si c bien un objet de type "Surv"
		if (classofY != "Surv") stop("Response must be a survival object")
	}

	ll <- attr(Terms, "term.labels")#liste des variables explicatives
#cluster(id) as.factor(dukes) as.factor(charlson) sex chemo terminal(death)

#=========================================================>

	mt <- attr(m, "terms") #m devient de class "formula" et "terms"

	X <- if (!is.empty.model(mt))model.matrix(mt, m, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables

#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
	strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
	cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
	num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
	vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()

	#booleen pour savoir si au moins une var depend du tps
	if (is.null(vartimedep)) timedep <- 0
	else timedep <- 1

	if (intcens & timedep) stop("You can not use time-varying effect covariates with interval censoring")
	if (intcens & cross.validation) stop("You can not do cross validation with interval censoring")
	if (intcens & logNormal) stop("It is currently impossible to fit a model with interval censoring and a log normal distribution of the frailties")
	if (timedep & logNormal) stop("You can not use time-varying effect covariates with a log normal distribution of the frailties")

	if(is.null(num.id)){
		joint.clust <- 1
	}else{
		joint.clust <- 0
		if (!joint) stop("num.id function can only be used with joint models")
	}
	
	subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

	if (length(subcluster)){
		ll <- ll[-grep("subcluster",ll)]
	}
	if (length(cluster)){
		ll_tmp <- ll[grep("cluster",ll)]
		ll <- ll[-grep("cluster",ll)]

		pos1 <- grep("r",unlist(strsplit(ll_tmp,split="")))[1]+2
		pos2 <- length(unlist(strsplit(ll_tmp,split="")))-1
		Names.cluster <- substr(ll_tmp,start=pos1,stop=pos2) # nom du cluster
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
		pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
		return(substr(x,start=pos1,stop=pos2))})
		occur <- rep(0,length(vec.factor))
	
		for(i in 1:length(vec.factor)){
			occur[i] = sum(vec.factor[i] == vect.fact)
		}
	}

#=========================================================>
	
	
	

	terminalEvent <- attr(Terms, "specials")$terminal #nbre de var qui sont en fonction de terminal()
	
	dropx <- NULL
	
	if (length(cluster) & Frailty == TRUE){
		tempc <- untangle.specials(Terms, "cluster", 1:10)
		ord <- attr(Terms, "order")[tempc$terms]
		if (any(ord > 1))stop("Cluster can not be used in an interaction")
		
		cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
		dropx <- tempc$terms
		uni.cluster<-unique(cluster)
	}else if (!length(cluster) & Frailty == TRUE){

		stop("grouping variable is needed")

	}else if (length(cluster) & Frailty == FALSE){
		stop("cluster not necessary for proportional hazard model")
	}
	else if (!length(cluster) & Frailty == FALSE){
		cluster <- 1:nrow(data) # valeurs inutiles pour un modele de Cox
		uni.cluster <- 1:nrow(data)
	}
	
	if (!missing(RandDist) & (Frailty == FALSE)){
		stop("RandDist not necessary for proportional hazard model")
	}

	if (length(num.id)){
		temppat <- untangle.specials(Terms, "num.id", 1:10)
		num.id <- m[,temppat$vars]
		dropx <- c(dropx,temppat$terms)
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
	
	#type <- attr(Y, "type")
	type <- typeofY

	if (type != "right" && type != "counting" && type != "interval" && type != "intervaltronc") { # Cox supporte desormais la censure par intervalle
		stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
	}

#	if ((type == "interval" || type == "interval2" || type == "intervaltronc") && intcens == FALSE) { # rajout
#		stop("You are trying to do interval censoring without intcens = TRUE")
#	}

	if (type != "counting" && recurrentAG) {
		stop("recurrentAG needs counting process formulation")
	}
	
	if (intcens == TRUE & recurrentAG == TRUE) {
		stop("recurrentAG invalid for interval censored data")
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
	Xlevels <- .getXlevels(newTerms, m)
	contr.save <- attr(X, 'contrasts')

	
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

	varnotdep <- colnames(X)[-grep("timedep",colnames(X))]
	vardep <- colnames(X)[grep("timedep",colnames(X))]
	vardep <- apply(matrix(vardep,ncol=1,nrow=length(vardep)),1,timedep.names)
	
	if (length(intersect(varnotdep,vardep)) != 0) {
		stop("A variable is both used as a constant and time-varying effect covariate")
	}

	nvartimedep <- length(vardep)

	filtretps <- rep(0,nvar)
	filtretps[grep("timedep",colnames(X))] <- 1

	var<-matrix(c(X),nrow=nrow(X),ncol=nvar) #matrix sans id et sans partie ex terminal(death)

	n<-nrow(X)


#add Alexandre 04/06/2012
#lire les donnees differemment si censure par intervalle
	if (intcens==TRUE) {
		if (type=="intervaltronc") {
			tt0 <- Y[,1]
			tt1 <- Y[,2]
			ttU <- Y[,3]
			cens <- Y[,4]
		} else {
			tt0 <- rep(0,n)
			tt1 <- Y[,1]
			ttU <- Y[,2]
			cens <- Y[,3]
			tt1[tt1==0] <- 0.1
		}
	} else {
		if (type=="right"){
			tt0 <- rep(0,n)
			tt1 <- Y[,1]
			cens <- Y[,2]
			ttU <- Y[,1] # rajouter quand meme dans frailPenal mais ne sera pas utilise
		} else {
			tt0 <- Y[,1]
			tt1 <- Y[,2]
			cens <- Y[,3]
			ttU <- Y[,2] # ne sera pas pris en compte dans le tri des temps de survie dans frailtypack.f90
		}                   # attention ne pas mettre de 0 sinon en cas de left trunc probleme dans la logV
	}

	if (min(cens)==0) cens.data<-1
	if (min(cens)==1 && max(cens)==1) cens.data<-0

	AG<-ifelse(recurrentAG,1,0)
	if (typeof == 0){
		crossVal<-ifelse(cross.validation,0,1)
	}


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
			indic.nb.int1 <- 0 # equals 0 for nb.int1 < 20
		}
		nbintervR <- nb.int1
		size1 <- 3*nbintervR
	}
	if ((typeof == 0) | (typeof == 2)) indic.nb.int1 <- 0

	if (sum(as.double(var))==0) nvar <- 0

	if (timedep==0){
		npbetatps <- 0
	}else{
		npbetatps <- (betaknots+betaorder-1)*nvartimedep
	}

	np <- switch(as.character(typeof),
		"0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + as.integer(Frailty)) + npbetatps,
		"1"=(as.integer(uni.strat) * nbintervR + nvar + as.integer(Frailty)) + npbetatps,
		"2"=(as.integer(uni.strat) * 2 + nvar + as.integer(Frailty)) + npbetatps)

	# traitement de l'initialisation du Beta rentre par l'utilisateur
	Beta <- rep(0,np)
	if (!missing(init.B)) {
		if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
		if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
		Beta <- c(rep(0,np-nvar),init.B)
	}
	if (!missing(init.Theta)) {
		if (!is.numeric(init.Theta)) stop("init.Theta must be numeric")
		if (Frailty==FALSE) stop("init.Theta does not exist in a Cox proportional hazard model")
		Beta[np-nvar] <- init.Theta
	}

		xSu1 <- rep(0,100)
		xSu2 <- rep(0,100)
		if (typeof==0){
			mt1 <- size1
		}else{
			mt1 <- 100
		}
		size2 <- mt1

		flush.console()
		if (print.times){
			ptm<-proc.time()
			cat("\n")
			cat("Be patient. The program is computing ... \n")
		}

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
				b=as.double(Beta),
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
				as.integer(intcens), # rajout
				as.double(ttU), # rajout
				logNormal=as.integer(logNormal),
				timedep=as.integer(timedep),
				as.integer(betaknots),
				as.integer(betaorder),
				as.integer(filtretps),
				BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep)),
				EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
				PACKAGE = "frailtypack") # 63
#AD:

    if (ans$istop == 4){
         warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
    }

    if (ans$istop == 2){
         warning("Model did not converge.")
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
#Al:
    fit$n.eventsbygrp <- table(cens,cluster)[2,]
    
    if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans[[24]]
    }else{
        fit$logLik <- ans[[24]]
    }

    if (Frailty) {
	fit$coef <- ans[[21]][(np - nvar - npbetatps + 1):np]
	if (logNormal == 0) fit$theta <- (ans[[21]][np - nvar - npbetatps])^2
	else fit$sigma2 <- (ans[[21]][np - nvar - npbetatps])^2
    }
    if (!Frailty) {
	if (logNormal == 0) fit$theta <- NULL
	else fit$sigma2 <- NULL
    }
    if (noVar1 == 1) {
      fit$coef <- NULL
    } 
    else
     {
	fit$coef <- ans[[21]][(np - nvar - npbetatps + 1):np]
	noms <- factor.names(colnames(X))
	if (timedep == 1){ # on enleve les parametres des B-splines qui ne serviront pas à l'utilisateur
		while (length(grep("timedep",noms))!=0){
			pos <- grep("timedep",noms)[1]
			noms <- noms[-pos]
			fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
		}
	}
	names(fit$coef) <- noms
     }

     temp1 <- matrix(ans[[22]], nrow = np, ncol = np)
     temp2 <- matrix(ans[[23]], nrow = np, ncol = np)
     if (Frailty) {
	fit$varTheta <- c(temp1[(np - nvar - npbetatps),(np - nvar - npbetatps)],temp2[(np - nvar - npbetatps),(np - nvar - npbetatps)])
	fit$varTheta <- ((2*ans[[21]][np - nvar - npbetatps])^2)*fit$varTheta # delta-method
     }
     
#AD:modification des dimensions des tableaux
    if(nvar > 0){

	 fit$varH <- temp1[(np - nvar - npbetatps + 1):np, (np - nvar - npbetatps + 1):np]
         fit$varHIH <- temp2[(np - nvar - npbetatps + 1):np, (np - nvar - npbetatps + 1):np]
	noms <- factor.names(colnames(X))
	if (timedep == 1){ # on enleve les variances des parametres des B-splines
		while (length(grep("timedep",noms))!=0){
			pos <- grep("timedep",noms)[1]
			noms <- noms[-pos]
			fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
			fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
		}
	}
   }
   
   fit$varHtotal <- temp1 # new Al: 20/06/13
   fit$varHIHtotal <- temp2
   

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

    fit$AG <- recurrentAG
    fit$intcens <- intcens # rajout
    fit$logNormal <- ans$logNormal

    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib
    fit$Names.data <- Names.data
    if (Frailty) fit$Names.cluster <- Names.cluster
    fit$Frailty <- Frailty
    if (Frailty){
	fit$martingale.res <- ans$martingale.res
	fit$frailty.pred <- ans$frailty.pred
	if (logNormal==0){
		fit$frailty.var <- ans$frailty.var
		fit$frailty.sd <- ans$frailty.sd
	}
    }else{
	fit$martingaleCox <- ans$martingaleCox
    }
    if (Frailty) fit$linear.pred <- ans$linear.pred[order(ordre)] # pour remettre dans le bon ordre
    else fit$linear.pred <- ans$linear.pred

  fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
  fit$nvartimedep <- nvartimedep

  fit$Names.vardep <- vardep
  
  fit$EPS <- ans$EPS

#
#========================= Test de Wald pour shared

    if(ans$istop==1){
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
	}

#===============================================
if (length(Xlevels) >0)fit$Xlevels <- Xlevels
    fit$contrasts <- contr.save
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
		if (classofY == "Surv")
		{
		if (!recurrentAG)
		{
			if(joint.clust==0){
				tempdc <- aggregate(tt1,by=list(num.id),FUN=sum)[,2]
				lignedc0 <- length(tempdc)
				tempdc <- cbind(rep(0,lignedc0),tempdc)
				clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				tt1.death <- 0
				tt0.death <- 0
			}else{
				tt1.death<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
				tt0.death<-rep(0,length(tt1.death))
				clusterdc <- 0
				lignedc0 <- 0
				tempdc <- 0
			}
		}else{
			if(joint.clust==0){
				#tempdc <- aggregate(tt1,by=list(num.id,cluster),FUN=sum)[,2]
				tempdc<-aggregate(tt1,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				lignedc0 <- length(tempdc)
				tempdc <- cbind(rep(0,lignedc0),tempdc)
				clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				tt1.death <- 0
				tt0.death <- 0

# 				tt0 <- aggregate(tt0,by=list(num.id),FUN=function(x) x[1])[,2]
# 				tt1 <- aggregate(tt1,by=list(num.id),FUN=function(x) x[1])[,2]
# 				cens <- aggregate(cens,by=list(num.id),FUN=function(x) x[1])[,2]
# 				cluster <- aggregate(cluster,by=list(num.id),FUN=function(x) x[1])[,2]
# 				if (!is.null(ncol(var))){ # si plus d'une variable explicative
# 					varAG<-aggregate(var[,1],by=list(num.id), FUN=function(x) x[1])[,2]
# 					if (ncol(var)>1){
# 						for (i in 2:ncol(var)){
# 							varAG.i<-aggregate(var[,i],by=list(num.id), FUN=function(x) x[1])[,2]
# 							varAG<-cbind(varAG,varAG.i)
# 						}
# 					}
# 					var<-varAG
# 				}else{
# 					var<-aggregate(var,by=list(num.id), FUN=function(x) x[1])[,2]
# 				}
# 				nobs <- n
# 				n <- length(tt0)

			}else{
				tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
				tt0.death<-rep(0,length(tt1.death))
				clusterdc <- 0
				lignedc0 <- 0
				tempdc <- 0

# 				tt0 <- aggregate(tt0,by=list(cluster),FUN=function(x) x[1])[,2]
# 				tt1 <- aggregate(tt1,by=list(cluster),FUN=function(x) x[1])[,2]
# 				ttU <- aggregate(ttU,by=list(cluster),FUN=function(x) x[1])[,2]
# 				cens <- aggregate(cens,by=list(cluster),FUN=function(x) x[1])[,2]
# 				if (!is.null(ncol(var))){ # si plus d'une variable explicative
# 					varAG<-aggregate(var[,1],by=list(cluster), FUN=function(x) x[1])[,2]
# 					if (ncol(var)>1){
# 						for (i in 2:ncol(var)){
# 							varAG.i<-aggregate(var[,i],by=list(cluster), FUN=function(x) x[1])[,2]
# 							varAG<-cbind(varAG,varAG.i)
# 						}
# 					}
# 					var<-varAG
# 				}else{
# 					var<-aggregate(var,by=list(cluster), FUN=function(x) x[1])[,2]
# 				}
# 				nobs <- n
# 				n <- length(tt0)

			}
		}
		}else{ # censure par intervalle
			if (recurrentAG == TRUE) stop("You can't fit joint models on interval-censored data with recurrentAG = TRUE")
			if(joint.clust==0){
				tempdc0 <- aggregate(tt0,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				tempdc <- aggregate(tt1,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				lignedc0 <- length(tempdc)
				#tempdc <- cbind(rep(0,lignedc0),tempdc)
				tempdc <- cbind(tempdc0,tempdc)
				clusterdc <- aggregate(cluster,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				tt1.death <- 0
				tt0.death <- 0

				# prendre en compte seulement un evenement pour le joint cluster
				tt0 <- aggregate(tt0,by=list(num.id),FUN=function(x) x[1])[,2]
				tt1 <- aggregate(tt1,by=list(num.id),FUN=function(x) x[1])[,2]
				ttU <- aggregate(ttU,by=list(num.id),FUN=function(x) x[1])[,2]
				cens <- aggregate(cens,by=list(num.id),FUN=function(x) x[1])[,2]
				cluster <- aggregate(cluster,by=list(num.id),FUN=function(x) x[1])[,2]
				if (!is.null(ncol(var))){ # si plus d'une variable explicative
				varAG<-aggregate(var[,1],by=list(num.id), FUN=function(x) x[1])[,2]
					if (ncol(var)>1){
						for (i in 2:ncol(var)){
							varAG.i<-aggregate(var[,i],by=list(num.id), FUN=function(x) x[1])[,2]
							varAG<-cbind(varAG,varAG.i)
						}
					}
					var<-varAG
				}else{
					var<-aggregate(var,by=list(num.id), FUN=function(x) x[1])[,2]
				}
				nobs <- n
				n <- length(tt0)

			}else{
				tt1.death<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
				tt0.death<-rep(0,length(tt1.death))
				clusterdc <- 0
				lignedc0 <- 0
				tempdc <- 0

				# prendre en compte seulement un evenement pour le joint
				tt0 <- aggregate(tt0,by=list(cluster),FUN=function(x) x[1])[,2]
				tt1 <- aggregate(tt1,by=list(cluster),FUN=function(x) x[1])[,2]
				ttU <- aggregate(ttU,by=list(cluster),FUN=function(x) x[1])[,2]
				cens <- aggregate(cens,by=list(cluster),FUN=function(x) x[1])[,2]
				if (!is.null(ncol(var))){ # si plus d'une variable explicative
					varAG<-aggregate(var[,1],by=list(cluster), FUN=function(x) x[1])[,2]
					if (ncol(var)>1){
						for (i in 2:ncol(var)){
							varAG.i<-aggregate(var[,i],by=list(cluster), FUN=function(x) x[1])[,2]
							varAG<-cbind(varAG,varAG.i)
						}
					}
					var<-varAG
				}else{
					var<-aggregate(var,by=list(cluster), FUN=function(x) x[1])[,2]
				}
				nobs <- n
				n <- length(tt0)
			}
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
			if(joint.clust==0){
				icdc00 <- aggregate(terminal,by=list(num.id),FUN=function(x) x[length(x)])[,2] #+aggregate(cens,by=list(num.id),FUN=function(x) x[length(x)])[,2]
				terminalEvent <- 0
			}else{
				terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
				icdc00 <- 0
			}
		}else{
			stop(" Joint frailty model miss specified ")
		}
#AD:


# terminalEvent might be 0-1
		if(joint.clust==0){
			if (!all(icdc00%in%c(1,0))){
				stop("terminal must contain a variable coded 0-1 and a non-factor variable")
			}
		}else{
			if (!all(terminalEvent%in%c(1,0))){
				stop("terminal must contain a variable coded 0-1 and a non-factor variable")
			}
		}
	
		m2 <- match.call(expand.dots = FALSE)

## AD: modified 20 06 2011, for no covariates on terminal event part
		if (missing(formula.terminalEvent)){
			m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <- m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m$... <- NULL
		}else{
			m2$formula.terminalEvent <- m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <- m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m$... <- NULL
		}

		
		m2$formula <- Terms2
		m2[[1]] <- as.name("model.frame")
		m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout
		

		match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA

		m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA
		
		if (!missing(formula.terminalEvent))newTerms2<-Terms2

#=========================================================>
		if (!missing(formula.terminalEvent)){
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
		}
	 
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
			Xlevels2 <- .getXlevels(newTerms2, m2)
			contr.save2 <- attr(X2, 'contrasts')
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

			nvar2 <- ncol(X2)

			vartimedep2 <- attr(Terms2, "specials")$timedep #nbre de var en fonction de timedep()

			# verifier qu'il y ait du timedep dans la deuxieme formule
			if (!is.null(vartimedep2)) timedep <- 1

			varnotdep2 <- colnames(X2)[-grep("timedep",colnames(X2))]
			vardep2 <- colnames(X2)[grep("timedep",colnames(X2))]
			vardep2 <- apply(matrix(vardep2,ncol=1,nrow=length(vardep2)),1,timedep.names)
			
			if (length(intersect(varnotdep2,vardep2)) != 0) {
				stop("A variable is both used as a constant and time-varying effect covariate in the formula of terminal event")
			}
		
			nvartimedep2 <- length(vardep2)

			filtretps2 <- rep(0,nvar2)
			filtretps2[grep("timedep",colnames(X2))] <- 1

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

			if(joint.clust==0){
				if (!is.null(ncol(vardc.temp))){
					vaxdc00<-aggregate(vardc.temp[,1],by=list(num.id), FUN=function(x) x[length(x)])[,2]  # num.id au lieu de cluster
					#print(head(vardc))
					
					if (ncol(vardc.temp)>1){
						for (i in 2:ncol(vardc.temp)){
							vaxdc00.i<-aggregate(vardc.temp[,i],by=list(num.id), FUN=function(x) x[length(x)])[,2]
							vaxdc00<-cbind(vaxdc00,vaxdc00.i)
						}
					}
				}else{
					vaxdc00<-aggregate(vardc.temp,by=list(num.id), FUN=function(x) x[length(x)])[,2]
				}
				vardc <- 0
			}else{
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
				vaxdc00 <- 0
			}
		}else{
			noVar2 <- 1
			vardc<-0
		}

	if ((classofY == "SurvIC") & (joint.clust==1) & (recurrentAG=FALSE)) { cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2] }
# 	cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2]
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
		if(joint.clust==1){
			if(is.null(nrow(vardc))){
				nvarEnd<-1
			}else{
				nvarEnd<-ncol(vardc)
			}
		}else{
			if(is.null(nrow(vaxdc00))){
				nvarEnd<-1
			}else{
				nvarEnd<-ncol(vaxdc00)
			}
		}
	}else{
		nvarEnd<-0
	}

	if (sum(as.double(var))==0) nvarRec <- 0
	if ((joint.clust==0) & sum(as.double(vaxdc00))==0) nvarEnd <- 0
	if ((joint.clust==1) & sum(as.double(vardc))==0) nvarEnd <- 0

	nvar<-nvarRec+nvarEnd

# ... end preparing data
#AD:
	effet <- 1
	indic_alpha <- 1
	if (!missing(Alpha)) { # new : joint more flexible alpha = 0
		if (Alpha=="none") indic_alpha <- 0
		else stop("Did you know that Alpha could only take 'none' as value ?")
	}
	
	nst <- 2
	
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

	if (timedep==0){
		npbetatps1 <- 0
		npbetatps2 <- 0
	}else{
		npbetatps1 <- (betaknots+betaorder-1)*nvartimedep
		npbetatps2 <- (betaknots+betaorder-1)*nvartimedep2
	}
	npbetatps <- npbetatps1 + npbetatps2

	np <- switch(as.character(typeof),
		"0"=(nst * (n.knots + 2) + nvarRec + nvarEnd + effet + indic_alpha + npbetatps),

		"1"=(nbintervR + nbintervDC + nvarRec + nvarEnd + effet + indic_alpha + npbetatps),

		"2"=(2*nst + nvarRec + nvarEnd + effet + indic_alpha + npbetatps))

	if (all(all.equal(as.numeric(cens),terminal)==T)){
		stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
	}

	# traitement de l'initialisation du Beta rentre par l'utilisateur
	Beta <- rep(0,np)
	if (!missing(init.B)) {
		if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
		if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
		Beta <- c(rep(0,np-nvar),init.B)
	}
	if (!missing(init.Theta)) {
		if (!is.numeric(init.Theta)) stop("init.Theta must be numeric")
		Beta[np-nvar-indic_alpha] <- init.Theta
	}
	if (!missing(init.Alpha)) {
		if (!missing(Alpha)) stop("You can not both initialize alpha parameter and fit a joint model without it")
		if (!is.numeric(init.Alpha)) stop("init.Alpha must be numeric")
		Beta[np-nvar] <- init.Alpha
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

	initialize <- 1

	npinit <- switch(as.character(typeof),
		"0"=((n.knots + 2) + nvarRec + effet),
		"1"=(nbintervR + nvarRec + nvarEnd + effet),
		"2"=(2 + nvarRec + effet))

	flush.console()
	if (print.times){
		ptm<-proc.time()
		cat("\n")
		cat("Be patient. The program is computing ... \n")
	}

	ans <- .Fortran("joint",
			as.integer(n),
			as.integer(length(uni.cluster)),
			as.integer(lignedc0),###
			as.integer(n.knots),
			k0=c(as.double(kappa1),as.double(kappa2)),
			as.double(tt0),
			as.double(tt1),
			as.integer(cens),
			as.integer(cluster),
			as.integer(clusterdc),###
			as.double(tt0.death),
			as.double(tt1.death),
			as.integer(terminalEvent),
			as.double(tempdc),###
			as.integer(icdc00),###
			as.integer(nvarRec),
			as.double(var),
			as.integer(nvarEnd),
			as.double(vardc),
			as.double(vaxdc00),###
			as.integer(noVar1),
			as.integer(noVar2),
			as.integer(recurrentAG),
			as.integer(maxit),
			np=as.integer(np),
			b=as.double(Beta),
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
			as.integer(typeof),
			as.integer(equidistant),
			as.integer(nbintervR),
			as.integer(nbintervDC),
			as.integer(c(size1,size2,mt11,mt12)),###
			ni=as.integer(0),
			cpt=as.integer(0),
			cpt.dc=as.integer(0),
			ier=as.integer(0),
			istop=as.integer(0),
			paraweib=as.double(rep(0,4)),
#			shape.weib=as.double(rep(0,2)),
#			scale.weib=as.double(rep(0,2)),
			MartinGale=as.double(matrix(0,nrow=length(uni.cluster),ncol=4)),###
			linear.pred=as.double(rep(0,n)),
			lineardc.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
			zi=as.double(rep(0,(n.knots+6))),
			time=as.double(rep(0,(nbintervR+1))),
			timedc=as.double(rep(0,(nbintervDC+1))),
# 			kendall=as.double(matrix(0,nrow=4,ncol=2)),
#			as.integer(initialize),
#			as.integer(npinit),
#			Bshared=as.double(rep(0,npinit)),
			linearpredG=as.double(rep(0,lignedc0)),
			joint.clust=as.integer(joint.clust),
			as.integer(c(intcens,indic_alpha)), # censure par intervalle, indic_alpha
			as.double(ttU),
			logNormal=as.integer(logNormal),
			paratps=as.integer(c(timedep,betaknots,betaorder)),
			as.integer(c(filtretps,filtretps2)),
			BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep)),
			BetaTpsMatDc=as.double(matrix(0,nrow=101,ncol=1+4*nvartimedep2)),
			EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
			PACKAGE = "frailtypack") # 65 = max


	MartinGale <- matrix(ans$MartinGale,nrow=as.integer(length(uni.cluster)),ncol=4)


    if (ans$istop == 4){
         warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
     }

    if (ans$istop == 2){
         warning("Model did not converge.")
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
    if (classofY == "SurvIC"){
        fit$n <- nobs
        if (typeofY == "intervaltronc") fit$indic.trunc <- 1
        else fit$indic.trunc <- 0
    }else{
        fit$n <- n
    }

    if (joint.clust == 0) fit$ind <- lignedc0
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

    if (logNormal == 0) fit$theta <- ans$b[np - nvar - npbetatps - indic_alpha]^2
    else fit$sigma2 <- ans$b[np - nvar - npbetatps - indic_alpha]^2

    if (indic_alpha == 1) fit$alpha <- ans$b[np - nvar - npbetatps]
    fit$npar <- np

#AD:
    if ((noVar1==1 & noVar2==1)) {
      fit$coef <- NULL
    }
    else
     {
       fit$coef <- ans$b[(np - nvar - npbetatps + 1):np]
	noms <- c(factor.names(colnames(X)),factor.names(colnames(X2)))
	if (timedep == 1){
		while (length(grep("timedep",noms))!=0){
			pos <- grep("timedep",noms)[1]
			noms <- noms[-pos]
			fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
		}
	}
	names(fit$coef) <- noms
#	if (missing(formula.terminalEvent)){
#	   names(fit$coef) <- c(factor.names(colnames(X)))
#	}else{
#          names(fit$coef) <- c(factor.names(colnames(X)), factor.names(colnames(X2)))
#      }
     }

#AD:
    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)
    
#Al:
    fit$varHIHtotal <- temp2
    
    fit$varH <- temp1[(np - nvar - npbetatps - indic_alpha):np, (np - nvar - npbetatps - indic_alpha):np]
    fit$varHIH <- temp2[(np - nvar - npbetatps - indic_alpha):np, (np - nvar - npbetatps - indic_alpha):np]
	if (indic_alpha == 1) noms <- c("theta","alpha",factor.names(colnames(X)),factor.names(colnames(X2)))
	else noms <- c("theta",factor.names(colnames(X)),factor.names(colnames(X2)))
	if (timedep == 1){ # on enleve les variances des parametres des B-splines
		while (length(grep("timedep",noms))!=0){
			pos <- grep("timedep",noms)[1]
			noms <- noms[-pos]
			fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
			fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
		}
	}
    fit$nvar<-c(nvarRec,nvarEnd)
    fit$nvarnotdep<-c(nvarRec-nvartimedep,nvarEnd-nvartimedep2)
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
    fit$shape.weib <- ans$paraweib[1:2]#ans$shape.weib
    fit$scale.weib <- ans$paraweib[3:4]#ans$scale.weib
#AD:

    # verif que les martingales ont ete bien calculees
    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
    if (Frailty){
	if (any(MartinGale[,1]==0)){
		fit$martingale.res <- msg
		fit$martingaledeath.res <- msg

		fit$frailty.pred <- msg
#		fit$frailty.var <- msg

		fit$linear.pred <- msg
		fit$lineardeath.pred <- msg
	}else{
		fit$martingale.res <- MartinGale[,1]#ans$martingale.res
		fit$martingaledeath.res <- MartinGale[,2]#ans$martingaledc.res
	
		fit$frailty.pred <- MartinGale[,3]#ans$frailty.pred
#		fit$frailty.var <- MartinGale[,4]#ans$frailty.var
	
		fit$linear.pred <- ans$linear.pred[order(ordre)]
		if (joint.clust==0){ fit$lineardeath.pred <- ans$linearpredG }
		else{ fit$lineardeath.pred <- ans$lineardc.pred }
	}
    }

#    if (joint.clust==0){
#        fit$kendall <- matrix(ans$kendall,nrow=4,ncol=2)
#    }
    fit$joint.clust <- ans$joint.clust
    fit$AG <- recurrentAG
    fit$intcens <- intcens # rajout
    fit$indic_alpha <- indic_alpha
    fit$logNormal <- ans$logNormal
    fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
    fit$BetaTpsMatDc <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedep2)
    fit$nvartimedep <- c(nvartimedep,nvartimedep2)

    fit$Names.vardep <- vardep
    fit$Names.vardepdc <- vardep2
    
    fit$EPS <- ans$EPS



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
	if (!missing(formula.terminalEvent)){
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
	}else{
		fit$global_chisq.test_d <- 0
	}
if (length(Xlevels) >0)	fit$Xlevels <- Xlevels
    fit$contrasts <- contr.save
if (length(Xlevels2) >0) fit$Xlevels2 <- Xlevels2
    fit$contrasts2 <- contr.save2
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
        if (logNormal == 1) stop("Nested model not implemented yet for log normal distribution of frailties")

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

	flush.console()
	if (print.times){
		ptm<-proc.time()
		cat("\n")
		cat("Be patient. The program is computing ... \n")
	}

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
		EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
		PACKAGE = "frailtypack") # 58

    if (ans$istop == 4){
	 warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
     }

    if (ans$istop == 2){
         warning("Model did not converge.")
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
       names(fit$coef) <- factor.names(colnames(X))
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
    fit$AG <- recurrentAG
    fit$EPS <- ans$EPS

 #   if (Frailty){

    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"
    if (any(ans$martingale.res==0)){

	fit$martingale.res <- msg
	fit$frailty.pred.group <- msg
	fit$frailty.pred.subgroup <- msg
	fit$linear.pred <- msg

    }else{

	fit$martingale.res <- ans$martingale.res
	fit$frailty.pred.group <- ans$frailty.pred.group

	nom1 <- paste("g_",c(1:ngg),sep="")
	nom2 <- paste("sub_g",c(1:maxng))

	frailty.pred.subgroup <- as.data.frame(matrix(round(ans$frailty.pred.subgroup,6),ncol=maxng))
	rownames(frailty.pred.subgroup) <- nom1
	colnames(frailty.pred.subgroup) <- nom2
	for (i in 1:ngg) {
		if (subgbyg[i] < max(subgbyg)) {
			frailty.pred.subgroup[i,(subgbyg[i]+1):max(subgbyg)] <- "."
		}
	}
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
	fit$linear.pred <- ans$linear.pred[order(ordre)]
    }
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
if (length(Xlevels) >0) fit$Xlevels <- Xlevels
    fit$contrasts <- contr.save
    attr(fit,"joint")<-joint
    attr(fit,"subcluster")<-TRUE
    class(fit) <- "nestedPenal"

 } # End NESTED MODEL

	if (print.times){
		cost<-proc.time()-ptm
		cat("The program took", round(cost[3],2), "seconds \n")
	}
 fit

}



