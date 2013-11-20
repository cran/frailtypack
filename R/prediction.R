
prediction <- function(fit, data, predTime, horizon, predMax, type="marginal", group, MC.sample=0){

	if (missing(fit)) stop("Need a fit")
	if ((class(fit)!="frailtyPenal") & (class(fit)!="jointPenal")) stop("The argument fit must be a frailtyPenal or jointPenal object")
	if (fit$istop != 1) stop("Attempting to do predictions with a wrong model")
	if (missing(data)) stop("Need data to do some predictions")
	if (missing(predTime) | missing(horizon) | missing(predMax)) stop("Need a first and a last time with an horizon to do predictions")
	if (horizon <= 0) stop("Horizon must be positive")
	if (predTime >= predMax) stop("Last time of predictions must be greater than first one")
	if (predTime < 0) stop("Be careful, negative time input")
	if ((MC.sample < 0) | (MC.sample > 1000))  stop("MC.sample needs to be positive integer up to 1000")
	
	if ((class(fit)=="jointPenal") & (!missing(type))) stop("No need for 'type' argument for predictions on a joint model")
	
	# seulement dans le cas du shared
	if (!(type %in% c("marginal","conditional"))) stop("Only 'marginal' or 'conditional' type of predictions can be specified")
	if (predMax > max(fit$x1)) stop("Prediction times cannot exceed maximum time of event")
	#if (!(predTime >= min(fit$x1))) stop("predtime must be in the right range")
	# mettre un warning quand une variable est un factor dans le fit et pas dans le datapred => source d'erreur
	
	if (MC.sample==0) ICproba <- FALSE
	else ICproba <- TRUE
	
	np <- fit$npar
	b <- fit$b
	typeof <- fit$typeof
	nva1 <- fit$nvarRec
	nva2 <- fit$nvarEnd
	ng <- fit$groups
	nst <- 2
	HIH <- fit$varHIHtotal

	# a definir meme si non utilise
	nz <- 1
	zi <- 0
	nbintervR <- 1
	nbintervDC <- 1
	time <- 0
	timedc <- 0
	
	if(typeof == 0){
		nz <- fit$n.knots.temp
		zi <- fit$zi
	}
	
	if(typeof == 1){
		nbintervR <- fit$nbintervR
		nbintervDC <- fit$nbintervDC
		time <- fit$time
		timedc <- fit$timedc
	}
	
	# nombre de predictions a faire pour chaque individu
	timeAll <- seq(predTime+horizon,predMax,by=horizon)
	ntimeAll <- length(timeAll)
	
	# recuperation des profils d'individus pour la prediction
	m <- fit$call
	m2 <- match.call()
	
	m$formula.terminalEvent <- m$Frailty <- m$joint <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$... <- NULL
	
	m[[1]] <- as.name("model.frame")
	m3 <- m # pour recuperer les donnees du dataset initial en plus
	m[[3]] <- as.name(m2$data)
	
	if (class(fit) == "jointPenal"){
		temp <- as.character(m$formula[[2]])
		if (temp[1]=="Surv"){
			if (length(temp) == 4) m$formula[[2]] <- as.name(temp[3])
			else if (length(temp) == 3) m$formula[[2]] <- as.name(temp[2])
			else stop("Wrong Surv function")
		}else{ # SurvIC
			if (length(temp) == 4) m$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],")"),collapse=" ")
			else if (length(temp) == 5) m$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],")"),collapse=" ")
			else stop("Wrong SurvIC function")
		}
		m$formula <- unlist(strsplit(deparse(m$formula)," "))
		m$formula <- gsub("\"","",m$formula)
		ter <- grep("terminal",m$formula)
		if (ter==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(ter,max(which(m$formula=="+")))],collapse=""))
		else m$formula <- as.formula(paste(m$formula[-ter],collapse=""))
		if (fit$joint.clust==0){
			m$formula <- unlist(strsplit(deparse(m$formula)," "))
			clus <- grep("cluster",m$formula)
			if (clus==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(clus,max(which(m$formula=="+")))],collapse=""))
			else m$formula <- as.formula(paste(m$formula[-clus],collapse=""))
		}
	}else{
		if (fit$Frailty){
			m$formula <- unlist(strsplit(deparse(m$formula)," "))
			clus <- grep("cluster",m$formula)
			if (clus==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(clus,max(which(m$formula=="+")))],collapse=""))
			else m$formula <- as.formula(paste(m$formula[-clus],collapse=""))
		}
		m$formula[[2]] <- NULL # pas besoin du Surv dans formula
	}
	
	dataset <- eval(m, sys.parent())
	dataset3 <- eval(m3, sys.parent())
	
	typeofY <- attr(model.extract(dataset3, "response"),"type")
	Y <- model.extract(dataset3, "response")
	
	if (typeofY=="right") tt1 <- Y[,1]
	else tt1 <- Y[,2]
	
	class(m$formula) <- "formula"
	special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")
	
	Terms <- terms(m$formula, special, data = data)
	
	m$formula <- Terms
	
	dropx <- NULL
	
	if (class(fit) == "jointPenal"){
		if (fit$joint.clust==1){ # joint classique
			tempc <- untangle.specials(Terms, "cluster", 1:10)
			dropx <- c(dropx,tempc$terms)
			cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
			uni.cluster <- unique(cluster)
			
			npred <- length(uni.cluster)
			nrec <- max(table(cluster))
			
			if (temp[1]=="Surv"){
				Y <- NULL
				for (i in uni.cluster) {
					temp <- model.extract(dataset, "response")[cluster==i]
					Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
				}
				predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
				trunctime <- rep(0,npred)
				lowertime <- rep(0,npred)
				uppertime <- rep(0,npred)
			}else{
				stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
				predtimerec <- matrix(0,nrow=npred)
				if (length(temp) == 4){ # pas troncature
					temp <- model.extract(dataset, "response")
					trunctime <- rep(0,npred)
					lowertime <- temp[,1]
					uppertime <- temp[,2]
				}
				if (length(temp) == 5){ # troncature
					temp <- model.extract(dataset, "response")
					trunctime <- temp[,1]
					lowertime <- temp[,2]
					uppertime <- temp[,3]
				}
			}
		}else{ # joint cluster
			tempnum <- untangle.specials(Terms, "num.id", 1:10)
			dropx <- c(dropx,tempnum$terms)
			num.id <- strata(dataset[, tempnum$vars], shortlabel = TRUE)
			uni.num.id <- unique(num.id)
			
			npred <- length(uni.num.id)
			nrec <- max(table(num.id))
			
			if (temp[1]=="Surv"){
				Y <- NULL
				for (i in uni.num.id) {
					temp <- model.extract(dataset, "response")[num.id==i]
					Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
				}
				predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
				trunctime <- rep(0,npred)
				lowertime <- rep(0,npred)
				uppertime <- rep(0,npred)
			}else{
				stop("Predictions not allowed for interval-censored yet...") # a enlever plus tard
				predtimerec <- matrix(0,nrow=npred)
				if (length(temp) == 4){ # pas troncature
					temp <- model.extract(dataset, "response")
					trunctime <- rep(0,npred)
					lowertime <- temp[,1]
					uppertime <- temp[,2]
				}
				if (length(temp) == 5){ # troncature
					temp <- model.extract(dataset, "response")
					trunctime <- temp[,1]
					lowertime <- temp[,2]
					uppertime <- temp[,3]
				}
			}
		}
	}else{
		if (fit$Frailty){
			class(m3$formula) <- "formula"
			Terms3 <- terms(m3$formula, special, data = data)
			m3$formula <- Terms3
			
			tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
			# je recupere le cluster du dataframe de depart et non pas des predictions
			cluster <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
			uni.cluster <- unique(cluster)
		}
	}
	
	if (!is.null(dropx)) newTerms <- Terms[-dropx]
	else newTerms <- Terms
	
	X <- model.matrix(newTerms, dataset)
	if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
	
	if (class(fit) == "jointPenal"){
	
		if (max(predtimerec) >= predTime) stop("Recurrences observed for each subject have to be less than predTime")
		
		if (fit$joint.clust==1) vaxpred <- aggregate(X,by=list(cluster),FUN=function(x) {x[1]})[,-1]
		else vaxpred <- aggregate(X,by=list(num.id),FUN=function(x) {x[1]})[,-1]
	
		# recuperation des variables partie deces
		m3 <- fit$call
		m2 <- match.call()
		
		m3$Frailty <- m3$joint <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa1 <- m3$kappa2 <- m3$maxit <- m3$hazard <- m3$nb.int1 <-m3$nb.int2 <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$... <- NULL
		
		m3$formula[[3]] <- m3$formula.terminalEvent[[2]]
		m3$formula.terminalEvent <- NULL
		m3[[1]] <- as.name("model.frame")
		m3[[3]] <- as.name(m2$data)
		
		temp <- as.character(m3$formula[[2]])
		if (temp[1]=="Surv"){
			if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
			else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
			else stop("Wrong Surv function")
		}else{ # SurvIC
			if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[2])
			else if (length(temp) == 5) m3$formula[[2]] <- as.name(temp[3])
			else stop("Wrong SurvIC function")
		}
		
		datasetdc <- eval(m3, sys.parent())
		
		class(m3$formula) <- "formula"
		special2 <- c("strata", "timedep")
		
		Terms2 <- terms(m3$formula, special2, data = data)
		
		X2 <- model.matrix(Terms2, datasetdc)
		if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
		if (fit$joint.clust==1) vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]
		else vaxdcpred <- aggregate(X2,by=list(num.id),FUN=function(x) {x[1]})[,-1]
		
		cat("\n")
		cat("Calculating the probabilities ... \n")
		
		ans <- .Fortran("predict",
				as.integer(np),
				as.double(b),
				as.integer(nz),
				as.integer(nbintervR),
				as.integer(nbintervDC),
				as.integer(nva1),
				as.integer(nva2),
				as.integer(nst),
				as.integer(typeof),
				as.double(zi),
				as.double(HIH),
				as.double(time),
				as.double(timedc),
				as.integer(ntimeAll),
				as.integer(npred),
				as.double(predTime),
				as.double(horizon),
				as.double(predtimerec),
				as.integer(nrec),
				as.double(as.matrix(vaxpred)),
				as.double(as.matrix(vaxdcpred)),
				proba1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				proba2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				proba3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probalow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probahigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probalow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probahigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probalow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				probahigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				icproba=as.integer(ICproba),
				as.integer(MC.sample),
				as.integer(fit$intcens),
				as.double(trunctime),
				as.double(lowertime),
				as.double(uppertime),
				PACKAGE = "frailtypack")
		
		out <- NULL
		out$call <- fit$call
		out$npred <- npred
		out$time <- timeAll
		if (fit$joint.clust==1) out$group <- uni.cluster
		else out$group <- uni.num.id
		if (!fit$intcens){
			out$proba1 <- matrix(ans$proba1,nrow=npred,ncol=ntimeAll)
			out$proba3 <- matrix(ans$proba3,nrow=npred,ncol=ntimeAll)
		}
		out$proba2 <- matrix(ans$proba2,nrow=npred,ncol=ntimeAll)
		out$icproba <- ICproba
		if (ICproba){
			if (!fit$intcens){
				out$probalow1 <- matrix(ans$probalow1,nrow=npred,ncol=ntimeAll)
				out$probahigh1 <- matrix(ans$probahigh1,nrow=npred,ncol=ntimeAll)
				out$probalow3 <- matrix(ans$probalow3,nrow=npred,ncol=ntimeAll)
				out$probahigh3 <- matrix(ans$probahigh3,nrow=npred,ncol=ntimeAll)
			}
			out$probalow2 <- matrix(ans$probalow2,nrow=npred,ncol=ntimeAll)
			out$probahigh2 <- matrix(ans$probahigh2,nrow=npred,ncol=ntimeAll)
		}
		out$joint.clust <- fit$joint.clust
		out$intcens <- fit$intcens
		
		cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")
		
		class(out) <- c("predJoint")
	}else{
		cat("\n")
		cat("Calculating the probabilities ... \n")
		
		expBX <- exp(X %*% fit$coef)
		
		# nombre de predictions a faire pour chaque individu
		# sequence <- seq(predTime,predMax,length=50)
		sequence <- seq(predTime+horizon,predMax,by=horizon)
		predMat <- NULL
		
		if (fit$Frailty){
			if (type=="marginal"){ ## marginal ##
				for (k in 1:nrow(data)){
					vect.survival.X <- survival(predTime,ObjFrailty=fit)**expBX[k] # sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
					vect.survival.X.horizon <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
					pred <- 1-((1+fit$theta*(-log(vect.survival.X)))/(1+fit$theta*(-log(vect.survival.X.horizon))))**(1/fit$theta)
					predMat <- cbind(predMat,pred)
				}
			}else{ ## conditional ##
				if (missing(group)) stop("For conditional predictions, a group need to be input")
				if (!(group %in% uni.cluster)) stop("Are you sure that the group is present in your cluster variable ?")
				
				for (k in 1:nrow(data)){
					vect.survival.X <- survival(predTime,ObjFrailty=fit)**expBX[k] # sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
					vect.survival.X.horizon <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
					pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[uni.cluster==group]
					predMat <- cbind(predMat,pred)
				}
			}
		}else{ ## for Cox model
			if (!missing(group)) stop("No need for a group to predict on a proportionnal hazard model")
			
			for (k in 1:nrow(data)){
				vect.survival.X <- survival(predTime,ObjFrailty=fit)**expBX[k] # sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
				vect.survival.X.horizon <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
				pred <- 1-(vect.survival.X.horizon/vect.survival.X)
				predMat <- cbind(predMat,pred)
			}
		}
		
		# -------------------------------------------------------- #
		# calcul des bornes de confiances (methode de Monte Carlo) #
		# -------------------------------------------------------- #
		
		if (ICproba){
			balea <- mvrnorm(MC.sample,fit$b,fit$varHtotal)
			#print(fit$b)
			#print(apply(balea,2,mean))
			#print(fit$b-apply(balea,2,mean))
			
			if (fit$Frailty) theta.mc <- balea[,fit$np-fit$nvar]^2
			
			aleaCoef <- balea[,(fit$np-fit$nvar+1):(fit$np)]
			expBX.mc <- exp(X %*% t(aleaCoef))
			
			# recuperation parametres de la fonction de risque/survie (splines,piecewise,weibull)
			if (fit$typeof == 0){
				para.mc <- balea[,1:(fit$n.knots+2)]^2
				if(fit$n.strat == 2) para.mc2 <- balea[,(fit$n.knots+3):(2*(fit$n.knots+2))]^2
				else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$n.knots+2)
			}else if (fit$typeof == 1){
				para.mc <- balea[,1:(fit$nbintervR)] # attention de ne pas éléever au carré
				if(fit$n.strat == 2) para.mc2 <- balea[,(fit$nbintervR+1):(2*fit$nbintervR)]
				else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$nbintervR)
			}else{
				para.mc <- balea[,1:2]^2
				if(fit$n.strat == 2) para.mc2 <- balea[,2:4]^2
				else para.mc2 <- matrix(0,nrow=MC.sample,ncol=2)
			}
			
			survival.mc <- function(t,ObjFrailty,para1,para2){ # dans les trois cas para1 et para2 seront traites differemment
				if (ObjFrailty$typeof == 0){ # splines
					nz <- ObjFrailty$n.knots
					zi <- ObjFrailty$zi
					res <- NULL
					nst <- ObjFrailty$n.strat
					out <- .Fortran("survival",as.double(t),as.double(para1),as.double(para2),as.integer(nz+2),
					as.double(zi),survival=as.double(c(0,0)),lam=as.double(c(0,0)),as.integer(nst),PACKAGE = "frailtypack") # lam ajoute suite aux modif de survival
					if(ObjFrailty$n.strat == 2){
						res <- c(res,out$survival)
					}else{
						res <- c(res,out$survival[1])
					}
					return(res)
				}
				if (ObjFrailty$typeof == 1){ # piecewise
					res <- NULL
					if (ObjFrailty$n.strat == 2) b <- c(para1,para2)
					else b <- para1
					time <- ObjFrailty$time
					out <- .Fortran("survival_cpm",as.double(t),as.double(b),
					as.integer(ObjFrailty$n.strat),as.integer(ObjFrailty$nbintervR),
					as.double(time),survival=as.double(c(0,0)),PACKAGE = "frailtypack")
					if(ObjFrailty$n.strat == 2){
						res <- c(res,out$survival)
					}else{
						res <- c(res,out$survival[1])
					}
					return(res)
				}
				if (ObjFrailty$typeof == 2){ # weibull
					res <- NULL
					sh1 <- para1[1]
					sc1 <- para1[2]
					res <- c(res,exp(-(t/sc1)^sh1))
					if(ObjFrailty$n.strat == 2){
						sh1 <- para2[1]
						sc1 <- para2[2]
						res <- c(res,exp(-(t/sc1)^sh1))
					}
					return(res)
				}
			}
			
			# calcul de la somme des risques cumules juste pour le groupe defini
			X3 <- model.matrix(newTerms, dataset3)
			if (ncol(X3) > 1) X3 <- X3[, -1, drop = FALSE]
			
			expBX3 <- exp(X3 %*% fit$coef)
			if ((fit$Frailty) & (type=="conditional")) res1 <- sum((-log(sapply(tt1[which(cluster==group)],survival,ObjFrailty=fit))) %*% expBX3[which(cluster==group)])
			
			predMatLow <- NULL
			predMatHigh <- NULL
			frailty.mc <- NULL

			if (fit$Frailty){
				if (type=="marginal"){ ## marginal ##
					for (k in 1:nrow(data)){
						realisations <- NULL
						for (i in 1:MC.sample){
							vect.survival.X <- survival.mc(predTime,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i] # sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							vect.survival.X.horizon <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							pred <- 1-((1+theta.mc[i]*(-log(vect.survival.X)))/(1+theta.mc[i]*(-log(vect.survival.X.horizon))))**(1/theta.mc[i])
							realisations <- cbind(realisations,pred)
						}
						predMatLow <- cbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
						predMatHigh <- cbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
					}
				}else{ ## conditional ##
					for (k in 1:nrow(data)){
						realisations <- NULL
						for (i in 1:MC.sample){
						
							if (k == 1) frailty.mc <- c(frailty.mc,rgamma(1,shape=fit$n.eventsbygrp[uni.cluster==group]+1/theta.mc[i],scale=1/(res1+1/theta.mc[i])))
							
							vect.survival.X <- survival.mc(predTime,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i] # sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							vect.survival.X.horizon <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							pred <- 1-(vect.survival.X.horizon/vect.survival.X)**frailty.mc[i]
							realisations <- cbind(realisations,pred)
						}
						predMatLow <- cbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
						predMatHigh <- cbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
					}
				}
			}else{ ## for a Cox model
				for (k in 1:nrow(data)){
					realisations <- NULL
					for (i in 1:MC.sample){
						vect.survival.X <- survival.mc(predTime,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i] # sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
						vect.survival.X.horizon <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
						pred <- 1-(vect.survival.X.horizon/vect.survival.X)
						realisations <- cbind(realisations,pred)
					}
					predMatLow <- cbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
					predMatHigh <- cbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
				}
			}
		}
		
		out <- NULL
		out$call <- fit$call
		out$time <- sequence
		out$pred <- predMat
		out$icproba <- ICproba
		if (ICproba){
			out$predLow <- predMatLow
			out$predHigh <- predMatHigh
		}
		if (fit$Frailty) out$type <- type
		out$horizon <- horizon
		if (type == "conditional") out$group <- group
		
		cat("Predictions done for",nrow(data),"subjects \n")
		
		class(out) <- "predFrailty"
	}
	out
}