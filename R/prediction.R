prediction <- function(fit, data, data.Longi, t, window, event = "Both", conditional=FALSE, MC.sample=0){
	# set.global <- function (x, value) { # Combinee a la fonction aggregate personnalisee, permet de tenir compte de variable dependantes du temps
		# x <- deparse(substitute(x))
		# assign(x, value, pos=.GlobalEnv)
	# }

	if (missing(fit)) stop("Need a fit")
	if ((class(fit)!="frailtyPenal") & (class(fit)!="jointPenal") & class(fit)!='longiPenal' & class(fit)!='trivPenal') stop("The argument fit must be a frailtyPenal or jointPenal object")
	if (fit$istop != 1) stop("Attempting to do predictions with a wrong model")
	
	if ((class(fit) == "jointPenal") && (fit$joint.clust == 0)) stop("Prediction method is not available for joint model for clustered data")

	if (missing(data)) stop("Need data to do some predictions")

	if (missing(data.Longi) & (class(fit)=="longiPenal") & (class(fit)=="trivPenal")) stop("Need data.Longi to do predictions")
	
	if (missing(t) | missing(window)) stop("Need times and a window to do predictions")
	if (length(t)!=1 & length(window)!=1) stop("t and window can not be vector both at the same time")
	if (is.unsorted(t)) stop("Last time of predictions must be greater than first one")
	if (any(t < 0)) stop("Be careful, negative time input")

	if ((class(fit) == "frailtyPenal") && ((fit$Frailty==TRUE) & (!missing (event)) && (event != "Recurrent"))){ # Prediction modele shared pour evenement repete
		stop("Only 'Recurrent' event is allowed for a shared frailty modeling of parameters")
	}

	if (any(window <= 0)) stop("Window must be positive")
	
	event.type <- charmatch(event, c("Both", "Terminal", "Recurrent"), nomatch = 0)
	if (event.type == 0) {
		stop("event must be 'Both', 'Terminal' or 'Recurrent'")
	}
	
	if ((MC.sample < 0) | (MC.sample > 1000))  stop("MC.sample needs to be positive integer up to 1000")

	if ((class(fit)=="jointPenal" | class(fit)=='longiPenal' | class(fit)=='trivPenal') & (conditional)) stop("No conditional prediction available on a joint model")

	if(class(fit)=='jointPenal' | class(fit)=='trivPenal'){
		if (max(t+window) > max(fit$xR)) stop("Prediction times cannot exceed maximum time of observation")
		if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")
	}
		
	if(class(fit)=='frailtyPenal'){
		if (max(t+window) > max(fit$x)) stop("Prediction times cannot exceed maximum time of observation")
	}

	if(class(fit)=='longiPenal'){
		if (max(t+window) > max(fit$xD)) stop("Prediction times cannot exceed maximum time of observation")
	}
	# seulement dans le cas du shared
	# if (missing(group)) type <- "marginal"
	# else type <- "conditional"

	# if (!(predTime >= min(fit$x1))) stop("predtime must be in the right range")
	# mettre un warning quand une variable est un factor dans le fit et pas dans le datapred => source d'erreur

	if (MC.sample==0) ICproba <- FALSE
	else ICproba <- TRUE

	np <- fit$npar
	b <- fit$b
	typeof <- fit$typeof
	nva1 <- fit$nvarRec
	nva2 <- fit$nvarEnd
	nva3 <- fit$nvarY
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
	moving.window <- FALSE
	if (length(t)==1) moving.window <- TRUE

	if (moving.window){
		predTime <- t
		timeAll <- t+window #seq(predTime+window,predMax,by=window)
		if (class(fit) == "jointPenal" | class(fit)== "trivPenal") window <- 0
	}else{
		predTime <- t[1]
		timeAll <- t+window
	}
	ntimeAll <- length(timeAll)
	formula_fit <- fit$formula

	# recuperation des profils d'individus pour la prediction
	m <- fit$call
	m2 <- match.call()

	m$formula.terminalEvent <- m$formula.LongitudinalData <- m$data.Longi <- m$random <- m$id  <- m$link <- m$left.censoring <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$init.Random <- m$init.Eta <- m$Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$...  <- NULL
		
	m[[1]] <- as.name("model.frame")
	m3 <- m # pour recuperer les donnees du dataset initial en plus
	m3$formula <- fit$formula
	m[[3]] <- as.name(m2$data)
			
	if (class(fit) == "jointPenal" | class(fit)=="trivPenal"){
		temp <- as.character(fit$formula[[2]])
		if (temp[1]=="Surv"){
			if (length(temp) == 4) fit$formula[[2]] <- paste(c("cbind(",temp[3],",",temp[4],")"),collapse=" ")
			else if (length(temp) == 3) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],")"),collapse=" ")
			else stop("Wrong Surv() function")
		}else{ # SurvIC
			if (length(temp) == 4) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],")"),collapse=" ")
			else if (length(temp) == 5) fit$formula[[2]] <- paste(c("cbind(",temp[2],",",temp[3],",",temp[4],",",temp[5],")"),collapse=" ")
			else stop("Wrong SurvIC() function")
		}
		fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
		fit$formula <- gsub("\"","",fit$formula)

		ter <- grep("terminal",fit$formula)
		if (ter==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(ter,max(which(fit$formula=="+")))],collapse=""))
		else m$formula <- as.formula(paste(fit$formula[-ter],collapse=""))

		if (fit$joint.clust==0){
			fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
			clus <- grep("cluster",fit$formula)
			if (clus==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(clus,max(which(fit$formula=="+")))],collapse=""))
			else m$formula <- as.formula(paste(fit$formula[-clus],collapse=""))
		}
	}else{
		if (fit$Frailty==TRUE ){
			if(event == 'Recurrent'){
				temp <- as.character(fit$formula[[2]])
				if (temp[1]=="Surv"){
					# if (length(temp) == 4) 
					fit$formula[[2]] <- paste(c("cbind(",temp[length(temp)-1],",",temp[length(temp)],")"),collapse=" ")
					# else if (length(temp) == 3) fit$formula[[2]] <- paste(c("cbind(",temp[length(temp)-1],",",temp[length(temp)],")"),collapse=" ")
					if (all(length(temp) != c(3,4))) stop("Wrong Surv() function")
				}
				fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
				fit$formula <- gsub("\"","",fit$formula)
				clus <- grep("cluster",fit$formula)
				m$formula <- as.formula(paste(fit$formula,collapse=""))
			}else{
				fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
				clus <- grep("cluster",fit$formula)
				# if (clus==length(fit$formula)) m$formula <- as.formula(paste(fit$formula[-c(clus,max(which(fit$formula=="+")))],collapse=""))
				# else m$formula <- as.formula(paste(fit$formula[-clus],collapse=""))				
				m$formula <- as.formula(paste(fit$formula,collapse=""))
			}
		}else if(class(fit)=='longiPenal'){
			fit$formula <- unlist(strsplit(deparse(fit$formula)," "))
			m$formula <- as.formula(paste(fit$formula,collapse=""))
		}else{#Cox
			m$formula <- fit$formula
		}
		if (!(fit$Frailty==TRUE && event == 'Recurrent')) m$formula[[2]] <- NULL # pas besoin du Surv dans formula, sauf si prediction pour evenement recurrent 
        # m$formula[[2]] <- NULL # pas besoin du Surv dans formula		
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
	fit$formula <- Terms

	dropx <- NULL	

	if (class(fit) == 'jointPenal' | class(fit) == 'trivPenal'){	
		if (fit$joint.clust==1){ # joint classique	
			tempc <- untangle.specials(Terms, "cluster", 1:10)		
			nameGrup <- substr(tempc$vars,9,nchar(tempc$vars)-1)
			
			dropx <- c(dropx,tempc$terms)
			cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
			uni.cluster <- unique(cluster)

			ic <- model.extract(dataset, "response")[,2]
			npred <- length(uni.cluster)
			nrec <- max(table(cluster[ic==1]))

			if (temp[1]=="Surv"){
				Y <- NULL
				for (i in uni.cluster) {
					temp <- model.extract(dataset, "response")[,1]
					temp <- temp[cluster==i & ic==1]
					Y <- c(Y,c(temp,rep(NA,nrec-length(temp))))
					#Y <- c(Y,c(temp,rep(0,nrec-length(temp))))
					# NA car permet de differencier les t=0 donnes par l utilisateur de ceux qui sont remplis par defaut
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

			ic <- model.extract(dataset, "response")[,2]
			npred <- length(uni.num.id)
			nrec <- max(table(num.id[ic==1]))
			
			if (temp[1]=="Surv"){
				Y <- NULL
				for (i in uni.num.id){
					temp <- model.extract(dataset, "response")[,1]
					temp <- temp[num.id==i & ic==1]
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
			#--Traitement donnees prediction
			tempc <- untangle.specials(Terms, "cluster", 1:10)
			dropx <- c(dropx,tempc$terms)
			cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)		
			uni.cluster <- unique(cluster)
			
			#--Traitement donnees fit$
			class(m3$formula) <- "formula"
			Terms3 <- terms(m3$formula, special, data = data)
			m3$formula <- Terms3
			tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
			clusterfit <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
			uni.clusterfit <- unique(clusterfit)
			
			#--Pour labelliser les lignes de l output
			nameGrup <- substr(tempc3$vars,9,nchar(tempc3$vars)-1)
			
			if(event == 'Recurrent'){
				# class(m3$formula) <- "formula"
				# Terms3 <- terms(m3$formula, special, data = data)
				# tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
				# clusterfit <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
				# uni.clusterfit <- unique(clusterfit)
				# nameGrup <- substr(tempc$vars,9,nchar(tempc$vars)-1)
				
				ic <- model.extract(dataset, "response")[,2]
				npred <- length(uni.cluster)
				nrec <- max(table(cluster[ic==1]))

				if (temp[1]=="Surv"){
					Y <- NULL
					for (i in uni.cluster) {
						temp <- model.extract(dataset, "response")[,1]
						temp <- temp[cluster==i & ic==1]
						Y <- c(Y,c(temp,rep(NA,nrec-length(temp))))
						# NA car permet de differencier les t=0 donnes par l utilisateur de ceux qui sont remplis par defaut
					}
					predtimerec <- matrix(Y,nrow=npred,byrow=TRUE)
                    predtimerectmp <- predtimerec 
                    predtimerectmp[which(is.na(predtimerectmp))] <- 0
				}
			} #else{
				# # class(m3$formula) <- "formula"
				# # Terms3 <- terms(m3$formula, special, data = data)
				# # m3$formula <- Terms3
				# # tempc3 <- untangle.specials(Terms3, "cluster", 1:10)
				# # cluster <- strata(dataset3[, tempc3$vars], shortlabel = TRUE)
				# # uni.cluster <- unique(cluster)
				# nameGrup <- substr(tempc3$vars,9,nchar(tempc3$vars)-1)
			# }
		}
	}
	if (!is.null(dropx)) newTerms <- Terms[-dropx]
	else newTerms <- Terms
	
	
	X <- model.matrix(newTerms, dataset)
	if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
	#####-----------------------------------------------------------------------------#####
	#####-&-&-&-&-&-&-&-&-&-&-&-Prediction for joint frailty model-&-&-&-&-&-&-&-&-&-&#####
	#####-----------------------------------------------------------------------------#####	
	if (class(fit) == "jointPenal"){	
		#--------ML 30-11-16...
		predtimerec <- predtimerec[order(unique(cluster)),]
		
		save(predtimerec, file = 'predtimerec.RData')
		save(predTime, file = 'predTime.RData')
		
		
		# listPrec <- NULL	
		# for (k in 1:nrow(predtimerec)){
			# tPrec <- which(predtimerec[k,] < predTime)
			# tPrec <- tPrec[length(tPrec)]
			# if (length(tPrec) == 0) tPrec <- 1 
			# listPrec <- c(listPrec,tPrec)
		# }
		
		taille = 0
		listPrec <- NULL  
		for (k in 1:nrow(predtimerec)){
			tPrec <- which(predtimerec[k,] < predTime)   
			if (length(tPrec) == 0) tPrec <- taille + 1 
			tPrec <- taille + length(tPrec)  
			
			rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
			if (length(rowTimes)==0) rowTimes <- 1
			taille = length(rowTimes)+taille
			listPrec <- c(listPrec,tPrec)                				
		}		
		
		predtimerec <- replace(predtimerec, is.na(predtimerec),0) #30-11-16	
		# listPrec <- rep(listPrec,ncol(X))		
		
		if (fit$joint.clust==1){#vaxpred <- aggregate(X,by=list(cluster),FUN=function(x) {x[1]})[,-1]		
			X <- X[order(cluster),]
			vaxpred <- X[listPrec,]			
		}else{ #vaxpred <- aggregate(X,by=list(num.id),FUN=function(x) {x[1]})[,-1]	
			X <- X[order(num.id),]
			vaxpred <- X[listPrec,]
		} 		
		#... ----ML 30-11-16		
		
		# recuperation des variables partie deces
		m3 <- fit$call
		m2 <- match.call()

		m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$... <- NULL

		m3$formula <- formula_fit
		m3$formula[[3]] <- fit$formula.terminalEvent[[2]]
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
		
		#------------ML: 01-12-16
		# listPrec <- rep(listPrec,ncol(X))
		if (fit$joint.clust==1){ #vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]			
			X2 <- X2[order(cluster),]
			vaxdcpred <- X2[listPrec,]
		}else { # vaxdcpred <- aggregate(X2,by=list(num.id),FUN=function(x) {x[1]})[,-1]
			X2 <- X2[order(num.id),]
			vaxdcpred <- X2[listPrec,]
		}	
		cat("\n")
		cat("Calculating the probabilities ... \n")
		#if(fit$logNormal==0){	#Myriam modifie le 18-08-16	
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
				as.integer(event.type),
				as.double(zi),
				as.double(HIH),
				as.double(time),
				as.double(timedc),
				as.integer(ntimeAll),
				as.integer(npred),
				as.double(predTime),
				as.double(window),
				as.double(predtimerec),
				as.integer(nrec),
				as.double(as.matrix(vaxpred)),
				as.double(as.matrix(vaxdcpred)),
				pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				pred1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				icproba=as.integer(ICproba),
				as.integer(MC.sample),
				as.integer(fit$intcens),
				as.double(trunctime),
				as.double(lowertime),
				as.double(uppertime),
				as.integer(moving.window),
				as.double(timeAll),
				as.integer(fit$logNormal),
				PACKAGE = "frailtypack") #43 arguments
				
		# Myriam 18-08-2016 Fusion des fichiers predict et predict_logN
		
		# }else{ #AK: joint log-normal
			# ans <- .Fortran("predict",
				# as.integer(np),
				# as.double(b),
				# as.integer(nz),
				# as.integer(nbintervR),
				# as.integer(nbintervDC),
				# as.integer(nva1),
				# as.integer(nva2),
				# as.integer(nst),
				# as.integer(typeof),
				# as.double(zi),
				# as.double(HIH),
				# as.double(time),
				# as.double(timedc),
				# as.integer(ntimeAll),
				# as.integer(npred),
				# as.double(predTime),
				# as.double(window),
				# #as.integer(event.type),
				# as.double(predtimerec),
				# as.integer(nrec),
				# as.double(as.matrix(vaxpred)),
				# as.double(as.matrix(vaxdcpred)),
				# pred1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# pred2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# pred3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# pred1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predlow1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predhigh1=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predlow2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predhigh2=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predlow3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predhigh3=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predlow1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# predhigh1_rec=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				# icproba=as.integer(ICproba),
				# as.integer(MC.sample),
				# as.integer(fit$intcens),
				# as.double(trunctime),
				# as.double(lowertime),
				# as.double(uppertime),
				# as.integer(moving.window),
				# as.double(timeAll), # 38
				# PACKAGE = "frailtypack")
		# }
		out <- NULL
		out$call <- match.call()
		out$name.fit <- match.call()[[2]]
		out$npred <- npred
		out$window <- window
		out$predtimerec <- predtimerec
		out$moving.window <- moving.window
		out$event <- event.type
		if (moving.window){
			out$x.time <- timeAll
			out$t <- predTime
		}else{
			out$x.time <- timeAll - window
		}
		if (fit$joint.clust==1) out$group <- uni.cluster[order(uni.cluster)]
		else out$group <- uni.num.id

		if (!fit$intcens){
			if ((event.type == 1) || (event.type == 2)){
				out$pred1 <- matrix(ans$pred1,nrow=npred,ncol=ntimeAll)
				rownames(out$pred1) <- paste(nameGrup,out$group)
				colnames(out$pred1) <- paste("time=", out$x.time)
				
				out$pred3 <- matrix(ans$pred3,nrow=npred,ncol=ntimeAll)
				rownames(out$pred3) <- paste(nameGrup,out$group)
				colnames(out$pred3) <- paste("time=", out$x.time)
			}
			
			if ((event.type == 1) || (event.type == 3)){
				out$pred1_rec <- matrix(ans$pred1_rec,nrow=npred,ncol=ntimeAll)
				rownames(out$pred1_rec) <- paste(nameGrup,out$group)
				colnames(out$pred1_rec) <- paste("time=", out$x.time)	
			}
		}
		if ((event.type == 1) || (event.type == 2)){
			out$pred2 <- matrix(ans$pred2,nrow=npred,ncol=ntimeAll)
			rownames(out$pred2) <- paste(nameGrup,out$group)
			colnames(out$pred2) <- paste("time=", out$x.time)
		}
		# Myriam : Modification de l'affichage des resultats (ajout des temps de prediction)
		out$icproba <- ICproba
		if (ICproba){
			if (!fit$intcens){ 
				if ((event.type == 1) || (event.type == 2)){
					out$predlow1 <- matrix(ans$predlow1,nrow=npred,ncol=ntimeAll)
					out$predhigh1 <- matrix(ans$predhigh1,nrow=npred,ncol=ntimeAll)
					rownames(out$predlow1) <- paste(nameGrup,out$group)
					colnames(out$predlow1) <- paste("time=", out$x.time)
					rownames(out$predhigh1) <- paste(nameGrup,out$group)
					colnames(out$predhigh1) <- paste("time=", out$x.time)
					
					out$predlow3 <- matrix(ans$predlow3,nrow=npred,ncol=ntimeAll)
					out$predhigh3 <- matrix(ans$predhigh3,nrow=npred,ncol=ntimeAll)
					rownames(out$predlow3) <- paste(nameGrup,out$group)
					colnames(out$predlow3) <- paste("time=", out$x.time)
					rownames(out$predhigh3) <- paste(nameGrup,out$group)
					colnames(out$predhigh3) <- paste("time=", out$x.time)
				}
				if ((event.type == 1) || (event.type == 3)){
					out$predlow1_rec <- matrix(ans$predlow1_rec,nrow=npred,ncol=ntimeAll)
					out$predhigh1_rec <- matrix(ans$predhigh1_rec,nrow=npred,ncol=ntimeAll)
					rownames(out$predlow1_rec) <- paste(nameGrup,out$group)
					colnames(out$predlow1_rec) <- paste("time=", out$x.time)
					rownames(out$predhigh1_rec) <- paste(nameGrup,out$group)
					colnames(out$predhigh1_rec) <- paste("time=", out$x.time)
				}
			}
			if ((event.type == 1) || (event.type == 2)){
				out$predlow2 <- matrix(ans$predlow2,nrow=npred,ncol=ntimeAll)
				out$predhigh2 <- matrix(ans$predhigh2,nrow=npred,ncol=ntimeAll)
				rownames(out$predlow2) <- paste(nameGrup,out$group)
				colnames(out$predlow2) <- paste("time=", out$x.time)
				rownames(out$predhigh2) <- paste(nameGrup,out$group)
				colnames(out$predhigh2) <- paste("time=", out$x.time)
			}
		}
		out$joint.clust <- fit$joint.clust
		out$intcens <- fit$intcens

		cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")
		class(out) <- c("predJoint")
	
	#####----------------------------------------------------------------------------#####
	#####-&-&-&-&-&-&-Prediction joint for longitudinal data and terminal event-&-&-&#####
	#####-&-&-&-&-&-&-		or longitudinal datas, recurrent events 		   -&-&-&#####
	#####-&-&-&-&-&-&-				and a terminal event					   -&-&-&#####
	#####----------------------------------------------------------------------------#####
	}else if(class(fit)=="longiPenal" | class(fit)=="trivPenal"){
		cat("\n")
		cat("Calculating the probabilities ... \n")
		
		if(class(fit)=="longiPenal"){	
			expBX <- exp(X %*% fit$coef[1:fit$nvarEnd])
		}else{
			#-----------ML:07-12-16
			taille = 0
			listPrec <- NULL  
			for (k in 1:nrow(predtimerec)){
				tPrec <- which(predtimerec[k,] < predTime)   
				if (length(tPrec) == 0) tPrec <- taille + 1 
				tPrec <- taille + length(tPrec) 
				
				rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
				if (length(rowTimes)==0) rowTimes <- 1
				taille = length(rowTimes)+taille
				listPrec <- c(listPrec,tPrec)
			}			
			predtimerec <- replace(predtimerec, is.na(predtimerec),0) #30-11-16
			listPrec <- listPrec[order(unique(cluster))]
		    vaxpred <- X[listPrec,]
			X <- X[order(cluster),]
			#-----------ML:07-12-16
									
			# recuperation des variables partie deces
			m3 <- fit$call
			m2 <- match.call()
			
			m3$formula.LongitudinalData <- m3$data.Longi <- m3$random <- m3$id <- m3$link <- m3$left.censoring <- m3$n.knots <- m3$recurrentAG <- m3$cross.validation <- m3$kappa <- m3$maxit <- m3$hazard <- m3$nb.int <- m3$RandDist <- m3$betaorder <- m3$betaknots <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$init.Theta <- m3$init.Alpha <- m3$Alpha <- m3$init.Random <- m3$init.Eta <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$... <- NULL
		
			m3$formula <- formula_fit
			m3$formula[[3]] <- fit$formula.terminalEvent[[2]]
			m3$formula.terminalEvent <- NULL
			m3[[1]] <- as.name("model.frame")
			m3[[3]] <- as.name(m2$data)

			temp <- as.character(m3$formula[[2]])

			if (length(temp) == 4) m3$formula[[2]] <- as.name(temp[3])
			else if (length(temp) == 3) m3$formula[[2]] <- as.name(temp[2])
			else stop("Wrong Surv function")
			   
			datasetdc <- eval(m3, sys.parent())
			class(m3$formula) <- "formula"
			special2 <- c("strata", "timedep")
			Terms2 <- terms(m3$formula, special2, data = data)
			X2 <- model.matrix(Terms2, datasetdc)
			if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]			

			vaxdcpred <- aggregate(X2,by=list(cluster),FUN=function(x) {x[1]})[,-1]			
		}
		# nombre de predictions a faire pour chaque individu
		if (moving.window){
			sequence2 <- t+window
			sequence <- rep(predTime,times=length(sequence2))
		}else{
			sequence <- t 
			sequence2 <- t+window
		}
		predMat <- NULL

		m2 <- fit$call
		m2$formula <-  m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard <- m2$nb.int <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$... <- NULL

		special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

		#========= Longitudinal Data preparation =========================
		class(fit$formula.LongitudinalData) <- "formula"
		TermsY <- terms(fit$formula.LongitudinalData, special, data = data.Longi)
		llY <- attr(TermsY, "term.labels")#liste des variables explicatives
		ord <- attr(TermsY, "order")
		
		#=========================================================>
		name.Y <- as.character(attr(TermsY, "variables")[[2]])	       
		
		if(class(fit)=="longiPenal") X <- X[order(unique(data.Longi$id)),]
		
		data.Longi <- data.Longi[order(data.Longi$id),]
		yy <- data.Longi[,which(names(data.Longi)==name.Y)]

		# on identifie les variables explicatives facteurs avec nombre de niveau plus que 2
		ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
		vec.factorY <- NULL
		vec.factorY <- c(vec.factorY,llY[ind.placeY])
		
		mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))

		# Fonction servant a prendre les termes entre "as.factor"
		vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
			if (length(grep("as.factor",x))>0){
				pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
				pos2 <- length(unlist(strsplit(x,split="")))-1
				return(substr(x,start=pos1,stop=pos2))
			}else{
				return(x)
			}
		})

		ind.placeY <- grep(paste(vec.factorY,collapse="|"),llY)
		if(is.factor(data.Longi[,names(data.Longi)==llY[1]])) X_L<- as.numeric(data.Longi[,names(data.Longi)==llY[1]])-1
		else X_L<- data.Longi[,names(data.Longi)==llY[1]]
		if(length(llY)>1){
			for(i in 2:length(llY)){
				if(is.factor(data.Longi[,names(data.Longi)==llY[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY[i]])-1)
				else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY[i]])
			}
		}
		#X_L<- data.Longi[,names(data.Longi)%in%(llY)]

		if(sum(ord)>length(ord)){
			for(i in 1:length(ord)){
				if(ord[i]>1){
					v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
					v2 <- strsplit(as.character(llY[i]),":")[[1]][2]
					if(is.factor(data.Longi[,names(data.Longi)==v1]) && length(levels(data.Longi[,names(data.Longi)==v1]))>2) stop("Interactions not allowed for factors with 3 or more levels (yet)")
					if(is.factor(data.Longi[,names(data.Longi)==v2]) && length(levels(data.Longi[,names(data.Longi)==v2]))>2) stop("Interactions not allowed for factors with 3 or more levels (yet)")
					if(is.factor(data.Longi[,names(data.Longi)==v1]) || !is.factor(data.Longi[,names(data.Longi)==v2])) {
						X_L <- cbind(X_L,(as.numeric(data.Longi[,names(data.Longi)==v1])-1)*data.Longi[,names(data.Longi)==v2])
						llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v1])[2],sep="")
					}else if (!is.factor(data.Longi[,names(data.Longi)==v1]) || is.factor(data.Longi[,names(data.Longi)==v2])) {
						X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*(as.numeric(data.Longi[,names(data.Longi)==v2])-1))
						llY[i]<-paste(llY[i],levels(data.Longi[,names(data.Longi)==v2])[2],sep="")
					}else{
						X_L <- cbind(X_L,data.Longi[,names(data.Longi)==v1]*data.Longi[,names(data.Longi)==v2])
					}
				}
			}
		}

		if(dim(X_L)[2]!=length(llY))stop("The variables in the longitudinal part must be in the data.Longi")
		X_L <- as.data.frame(X_L)
		names(X_L) <- llY
		Intercept <- rep(1,dim(X_L)[1])
		if(fit$intercept)X_L <- cbind(Intercept,X_L)
		X_Lall<- X_L
		"%+%"<- function(x,y) paste(x,y,sep="")
		if(length(vec.factorY) > 0){
			for(i in 1:length(vec.factorY)){
				X_L <- cbind(X_L[,-(which(names(X_L)==vec.factorY[i]))],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
			}
			vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]
			occurY <- rep(0,length(vec.factorY))
			for(i in 1:length(vec.factorY)){
				#occur[i] <- sum(vec.factor[i] == vect.fact)
				occurY[i] <- length(grep(vec.factorY[i],vect.factY))
			}
		}

		if (ncol(X_L) == 0){
			noVarY <- 1
		}else{
			noVarY <- 0
		}
		#=========================================================>
		clusterY <- data.Longi$id
		maxy_rep <- max(table(clusterY))
		uni.cluster<-as.factor(unique(clusterY))
		npred <- length(uni.cluster)
		nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:
		varY <- as.matrix(sapply(X_L, as.numeric))
		
		#=======================================>
		#======= Construction du vecteur des indicatrice
		if(length(vec.factorY) > 0){
			#ind.place <- ind.place -1	
			k <- 0
			for(i in 1:length(vec.factorY)){
				ind.placeY[i] <- ind.placeY[i]+k
				k <- k + occurY[i]-1		
			}
		}
		if(fit$link=="Random-effects")link <- 1
		if(fit$link=="Current-level") link <- 2
		if(fit$leftCensoring==FALSE){
			s_cag_id = 0
			s_cag = 0
		}else{
			s_cag_id = 1
			s_cag = fit$leftCensoring.threshold
		}				
		
		if(class(fit)=="longiPenal"){
			ans <- .Fortran("predict_biv",
				as.integer(np),
				as.double(b),
				as.integer(nz),
				as.integer(nva2),
				as.integer(nva3),
				as.integer(fit$ne_re),
				as.integer(fit$netadc),
				as.integer(link),
				as.integer(nst),
				as.integer(typeof),
				as.double(zi),
				as.double(HIH),
				as.integer(ntimeAll),
				as.integer(npred),
				as.double(predTime),
				as.double(window),
				as.integer(fit$max_rep),
				as.double(yy),
				as.double(as.matrix(X)),
				as.double(as.matrix(varY)),
				as.integer(clusterY),
				as.integer(unique(clusterY)),
				as.integer(length(clusterY)),
				pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				icproba=as.integer(ICproba),
				as.integer(MC.sample),
				as.integer(moving.window),
				as.double(timeAll),
				as.integer(s_cag_id),
				as.double(s_cag),
				PACKAGE = "frailtypack")#32 arguments

			predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)
			predMatLow <- matrix(ans$predlow,nrow=nrow(data),ncol=ntimeAll)
			predMatHigh <- matrix(ans$predhigh,nrow=nrow(data),ncol=ntimeAll)

			out <- NULL
			out$call <- match.call()
			out$name.fit <- match.call()[[2]]
			out$npred <- npred
			out$moving.window <- moving.window
			if (moving.window){
				out$x.time <- sequence2
				out$t <- predTime
			}else{
				out$x.time <- sequence
			}
			out$group <- uni.cluster
			out$pred <- predMat
			colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))
			#rownames(out$pred) <- paste("ind",1:out$npred)
			rownames(out$pred) <- paste("ind",unique(clusterY))

			out$icproba <- ICproba
			if (ICproba){
				out$predLow <- predMatLow
				out$predHigh <- predMatHigh
				colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
				# rownames(out$predLow) <- paste("ind",1:out$npred)
				rownames(out$predLow) <- paste("ind",unique(clusterY))
				colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
				# rownames(out$predHigh) <- paste("ind",1:out$npred)
				rownames(out$predHigh) <- paste("ind",unique(clusterY))
			}
			out$window <- window
			out$trivariate <- FALSE	
			
		}else if(class(fit)=="trivPenal"){             
			predtimerec <- predtimerec[order(unique(cluster)),]	
						
			ans <- .Fortran("predict_tri",
				as.integer(np),
				as.double(b),
				as.integer(nz),
				as.integer(nva1),
				as.integer(nva2),
				as.integer(nva3),
				as.integer(fit$ne_re),
				as.integer(fit$netar),
				as.integer(fit$netadc),
				as.integer(link),
				as.integer(nst),
				as.integer(typeof),
				as.double(zi),
				as.double(HIH),
				as.integer(ntimeAll),
				as.integer(npred),
				as.double(predTime),
				as.double(window),
				as.double(predtimerec),
				as.integer(nrec),
				as.integer(fit$max_rep),
				as.double(yy),
				as.double(as.matrix(vaxpred)),
				as.double(as.matrix(X)),
				as.double(as.matrix(varY)),
				as.integer(clusterY),
				as.integer(unique(clusterY)),
				as.integer(length(clusterY)),
				as.integer(npred),
				pred=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predlow=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				predhigh=as.double(matrix(0,nrow=npred,ncol=ntimeAll)),
				icproba=as.integer(ICproba),
				as.integer(MC.sample),
				as.integer(moving.window),
				as.double(timeAll),
				as.integer(s_cag_id),
				as.double(s_cag),
				PACKAGE = "frailtypack") #38 arguments
				
			out <- NULL
			out$call <- match.call()
			out$name.fit <- match.call()[[2]]
			out$npred <- npred
			out$window <- window
			out$predtimerec <- predtimerec
			out$moving.window <- moving.window
			if (moving.window){
				out$x.time <- timeAll
				out$t <- predTime
			}else{
				out$x.time <- timeAll - window
			} 
			out$group <- uni.cluster
			out$pred <- matrix(ans$pred,nrow=npred,ncol=ntimeAll)
			rownames(out$pred) <- paste("ind",out$group)
			colnames(out$pred) <- c("times",rep(" ",ntimeAll-1))
			
			out$icproba <- ICproba
			if (ICproba){
				out$predLow <- matrix(ans$predlow,nrow=npred,ncol=ntimeAll)
				out$predHigh <- matrix(ans$predhigh,nrow=npred,ncol=ntimeAll)
				rownames(out$predLow) <- paste("ind",out$group)
				colnames(out$predLow) <- c("times",rep(" ",ntimeAll-1))
				rownames(out$predHigh) <- paste("ind",out$group)
				colnames(out$predHigh) <- c("times",rep(" ",ntimeAll-1))
			}
			out$trivariate <- TRUE
		}
		cat("Predictions done for",npred,"subjects and",ntimeAll,"times \n")
		class(out) <- "predLongi"
		
	####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
	####-*-*-*-*-Prediction pour un modele Shared -*-*-*-*-####
	####*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*####
	}else if(class(fit)=="frailtyPenal"){	
		cat("\n")
		cat("Calculating the probabilities ... \n")		
		
		if (event == 'Recurrent'){	
			taille = 0
			listPrec <- NULL  
			for (k in 1:nrow(predtimerec)){
				tPrec <- which(predtimerec[k,] < predTime)   
				if (length(tPrec) == 0) tPrec <- taille + 1 
				tPrec <- taille + length(tPrec)    
				
				rowTimes <- predtimerec[k,][which(!is.na(predtimerec[k,]))]
				if (length(rowTimes)==0) rowTimes <- 1
				taille = length(rowTimes)+taille
				listPrec <- c(listPrec,tPrec)                				
			}			
			# X <- as.matrix(aggregate(X,by=list(cluster),FUN=function(x){x[1]})[,-1])
			X <- X[listPrec,]
			if (length(unique(cluster)) > 1) X <- X[order(unique(cluster)),]
		}
		expBX <- exp(X %*% fit$coef)	

		# nombre de predictions a faire pour chaque individu
		if (moving.window){ # 2 facons differentes de faire des predictions, soit h evolue, soit t evolue
			sequence2 <- t+window #seq(predTime+window,predMax,by=window)
			sequence <- rep(predTime,times=length(sequence2))
		}else{
			sequence <- t #seq(predTime,predMax,length=50)
			sequence2 <- t+window #sequence+window
		}
		predMat <- NULL

		if (fit$Frailty){
			###############################
			###   Prediction marginale  ###
			###############################
			if (!conditional){	
				#======================================================#
				#************ML: Prediction pour evenement recurrent****
				#======================================================#
				if(event == 'Recurrent'){				
					mat.survival.X <- NULL
					mat.survival.X.horizon <- NULL
					mat.survival.LastRec <- NULL
					npred0 <- nrow(predtimerec) #nb subjects
					nbrec <- rep(0,npred0)
					
					for (k in 1:npred0){					
						# vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
						vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)
						# vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
						vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)
						mat.survival.X <- rbind(mat.survival.X,vect.survival.X)
						mat.survival.X.horizon <- rbind(mat.survival.X.horizon,vect.survival.X.horizon)					
						
						recurr <- which(!is.na(predtimerec[k,which(predtimerec[k,] <= predTime)]))
						nbrec[k] <- length(recurr)						
												
						if(length(recurr) == 0) LastRec <- 0
						else if (length(recurr) > 1) LastRec <- predtimerec[k, recurr[-c(1:length(recurr)-1)]]
						else LastRec <- predtimerec[k,recurr]
						
						if((length(recurr) == 1)&&(LastRec == 0)) nbrec[k] <- 0
						
						# vect.survival.LastRec <- survival(LastRec,fit)**expBX[k]
						if (LastRec == 0) vect.survival.LastRec <- 1
						else vect.survival.LastRec <- survival(LastRec,fit)
						mat.survival.LastRec <- rbind(mat.survival.LastRec,vect.survival.LastRec)
					}					
					######### Distribution gamma #########			
					if(fit$logNormal==0) variance <- fit$theta
					######### Distribution LogNormale #########
					else variance <- fit$sigma2	
					
					nbrec <- nbrec[order(uni.cluster)]
					mat.survival.LastRec[,1] <- mat.survival.LastRec[order(uni.cluster)]	
					
					ans <- .Fortran("predict_Recurr_Sha",
						as.integer(fit$logNormal),
						as.integer(npred0),
						as.double(mat.survival.X),
						as.double(mat.survival.X.horizon),
						as.double(mat.survival.LastRec),
						as.double(expBX),
						as.double(variance),
						pred=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
						as.integer(nbrec),
						as.integer(ntimeAll),
						as.integer(0),						
						as.integer(MC.sample),
						as.double(rep(0,MC.sample)),
						as.double(matrix(0,nrow=npred0*MC.sample,ncol=ntimeAll)),
						as.double(matrix(0,nrow=npred0*MC.sample,ncol=ntimeAll)),
						as.double(rep(0,nrow=npred0*MC.sample)),
						as.double(rep(0,nrow=npred0*MC.sample)),
						predlow1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
						predhigh1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
						PACKAGE = "frailtypack")#19 arguments
 						
					predMat <- matrix(ans$pred,nrow=npred0,ncol=ntimeAll)
					
				#=================================================#
				#*         Prediction pour donnees groupees       *
				#=================================================#
				}else{				    
					######### Distribution gamma #########			
					if(fit$logNormal==0){	
						for (k in 1:nrow(data)){
							vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
							vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
							pred <- 1-((1+fit$theta*(-log(vect.survival.X)))/(1+fit$theta*(-log(vect.survival.X.horizon))))**(1/fit$theta)
							predMat <- rbind(predMat,pred)
						}			
					}else{
						######### AK: Distribution LogNormale #########				
						mat.survival.X <- NULL
						mat.survival.X.horizon <- NULL
						for (k in 1:nrow(data)){						
							#vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k] #ML : Survie doit tenir compte du terme de fragilite 
							vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)
							#vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]	#ML
							vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)
							mat.survival.X <- rbind(mat.survival.X,vect.survival.X)
							mat.survival.X.horizon <- rbind(mat.survival.X.horizon,vect.survival.X.horizon)							
						}								
						ans <- .Fortran("predict_LogN_sha",
							as.integer(nrow(data)),
							as.double(mat.survival.X),
							as.double(mat.survival.X.horizon),
							as.double(expBX),
							as.double(fit$sigma2),
							pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
							as.integer(0),
							as.integer(ntimeAll),
							as.integer(MC.sample),
							as.double(rep(0,MC.sample)),
							as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
							as.double(matrix(0,nrow=nrow(data)*MC.sample,ncol=ntimeAll)),
							predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
							predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
							PACKAGE = "frailtypack")	#14 arguments							
						predMat <- matrix(ans$pred,nrow=nrow(data),ncol=ntimeAll)
					}
				}				
			###############################
			###Prediction conditionnelle###
			###############################
			}else{ 			
				if (event =="Recurrent"){
					if (is.null(nrow(X))){						
						if (!(unique(cluster) %in% uni.clusterfit)) stop("Are you sure that the group",unique(cluster)," is present in your cluster variable ?")			
						vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX
						vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX
						pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[uni.clusterfit==as.integer(unique(cluster))]
						predMat <- rbind(predMat,pred)
					}else{
						uni.cluster <- uni.cluster[order(uni.cluster)]
						for (k in 1:nrow(X)){
							# if (!(group %in% uni.clusterfit)) stop("Are you sure that the group",group," is present in your cluster variable ?")					
							if (!(uni.cluster[k] %in% uni.clusterfit)) stop("Are you sure that the group",uni.cluster[k]," is present in your cluster variable ?")			
							vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
							vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]
							pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[uni.clusterfit==as.integer(uni.cluster[k])]
							predMat <- rbind(predMat,pred)
						}
					}
				}else{
					cluster <- as.integer(as.vector(cluster))														
					if (any(!(uni.cluster %in% uni.clusterfit))) stop("Are you sure that the group", uni.cluster, "is present in your cluster variable ?")					
					for (i in 1:nrow(data)){										
						vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[i]
						vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[i]
						if (fit$logNormal==0){ # Gamma distribution
							pred <- 1-(vect.survival.X.horizon/vect.survival.X)**fit$frailty.pred[uni.clusterfit==cluster[i]]
						}else{ #AK: Normal distribution
							pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(fit$frailty.pred[uni.clusterfit==cluster[i]])
						}
						predMat <- rbind(predMat,pred)
					}
				}
			}
	    #*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
		#-*-*-*-*-Pour un modele de Cox-*-*-*-*-#
		#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#		
		}else{
			# if (!missing(group)) stop("No need for a group to predict on a proportionnal hazard model")
			for (k in 1:nrow(data)){
				vect.survival.X <- sapply(sequence,FUN=survival,ObjFrailty=fit)**expBX[k]
				vect.survival.X.horizon <- sapply(sequence2,FUN=survival,ObjFrailty=fit)**expBX[k]				
				pred <- 1-(vect.survival.X.horizon/vect.survival.X)
				predMat <- rbind(predMat,pred)
			}
		}		
		# -------------------------------------------------------- #
		# calcul des bornes de confiances (methode de Monte Carlo) #
		# -------------------------------------------------------- #
		if (ICproba){
			balea <- mvrnorm(MC.sample,fit$b,fit$varHtotal)
			if (fit$Frailty){ #AK: For Gamma we have variance theta and for Normal we have variance sigma2
				if(fit$logNormal==0)theta.mc <- balea[,fit$np-fit$nvar]^2
				if(fit$logNormal==1)sigma2.mc <- balea[,fit$np-fit$nvar]^2
			}
			aleaCoef <- balea[,(fit$np-fit$nvar+1):(fit$np)]		
			expBX.mc <- exp(X %*% t(aleaCoef))
					
			# recuperation parametres de la fonction de risque/survie (splines,piecewise,weibull)
			if (fit$typeof == 0){ #Splines
				para.mc <- balea[,1:(fit$n.knots+2)]^2
				if(fit$n.strat == 2) para.mc2 <- balea[,(fit$n.knots+3):(2*(fit$n.knots+2))]^2
				else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$n.knots+2)				
			}else if (fit$typeof == 1){ #Piecewise
				para.mc <- balea[,1:(fit$nbintervR)] # attention de ne pas elever au carre
				if(fit$n.strat == 2) para.mc2 <- balea[,(fit$nbintervR+1):(2*fit$nbintervR)]
				else para.mc2 <- matrix(0,nrow=MC.sample,ncol=fit$nbintervR)				
			}else{ #Weibull
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
					
					out <- .Fortran("survival",
						as.double(t),
						as.double(para1),
						as.double(para2),
						as.integer(nz+2),
						as.double(zi),
						survival=as.double(c(0,0)),
						lam=as.double(c(0,0)),
						as.integer(nst),# lam ajoute suite aux modif de survival
						PACKAGE = "frailtypack") #8 arguments
						
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
					
					out <- .Fortran("survival_cpm",
						as.double(t),
						as.double(b),
						as.integer(ObjFrailty$n.strat),
						as.integer(ObjFrailty$nbintervR),
						as.double(time),
						survival=as.double(c(0,0)),
						PACKAGE = "frailtypack"
					) #6 arguments
					
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
				
			}# end of survival.mc function
			
			# calcul de la somme des risques cumules juste pour le groupe defini
			X3 <- model.matrix(newTerms, dataset3)
			if (ncol(X3) > 1) X3 <- X3[, -1, drop = FALSE]

			expBX3 <- exp(X3 %*% fit$coef)
			
			res1 <- sum((-log(sapply(tt1[which(clusterfit==5)],survival,ObjFrailty=fit))) %*% expBX3[which(clusterfit==5)])		
		
			predMatLow <- NULL
			predMatHigh <- NULL
			frailty.mc <- NULL			
			#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
			#-*-*-*-*-Pour un modele Shared-*-*-*-*-#
			#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
			if (fit$Frailty){
				###############################
				###   Prediction marginale  ###
				###############################
				if (!conditional){   
					#======================================================#
					#************ML: Prediction pour evenement recurrent****
					#======================================================#
					if(event == 'Recurrent'){	
						mat.survival.X.mc <- NULL
						mat.survival.X.horizon.mc <- NULL
						mat.survival.LastRec.mc <- NULL	
													
						for(i in 1:MC.sample){
							mat.survival.X.samp <- NULL
							mat.survival.X.horizon.samp <- NULL
							mat.survival.LastRec.samp <- NULL							
							for(k in 1:npred0){
								# vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
								vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
								# vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
								vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
								mat.survival.X.samp <- rbind(mat.survival.X.samp,vect.survival.X.samp)
								mat.survival.X.horizon.samp <- rbind(mat.survival.X.horizon.samp,vect.survival.X.horizon.samp)								
								
								recurr <- which(!is.na(predtimerec[k,which(predtimerec[k,] <= predTime)]))
								nbrec[k] <- length(recurr)
																
								if(length(recurr) == 0) LastRec <- 0
								else if (length(recurr) > 1) LastRec <- predtimerec[k,recurr[-1]]
								else LastRec <- predtimerec[k,recurr]	
								# vect.survival.LastRec.samp <- sapply(LastRec,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
								vect.survival.LastRec.samp <- sapply(LastRec,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
								mat.survival.LastRec.samp <- cbind(mat.survival.LastRec.samp,vect.survival.LastRec.samp)									
							}									
							mat.survival.X.mc <- rbind(mat.survival.X.mc,mat.survival.X.samp)
							mat.survival.X.horizon.mc <- rbind(mat.survival.X.horizon.mc,mat.survival.X.horizon.samp)
							mat.survival.LastRec.mc <- rbind(mat.survival.LastRec.mc,mat.survival.LastRec.samp)
						}
												
						######### Distribution gamma #########				
						if(fit$logNormal==0){ 
							variance <- fit$theta
							variance.mc <- theta.mc
						######### Distribution LogNormale #########
						}else{
							variance <- fit$sigma2
							variance.mc <- sigma2.mc
						}	
						
						ans <- .Fortran("predict_Recurr_Sha",
							as.integer(fit$logNormal),
							as.integer(npred0),
							as.double(mat.survival.X),
							as.double(mat.survival.X.horizon),
							as.double(mat.survival.LastRec),
							as.double(expBX),
							as.double(variance),
							pred=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
							as.integer(nbrec),
							as.integer(ntimeAll),
							as.integer(1),						
							as.integer(MC.sample),
							as.double(variance.mc),
							as.double(mat.survival.X.mc),
							as.double(mat.survival.X.horizon.mc),
							as.double(mat.survival.LastRec.mc),
							as.double(expBX.mc),
							predlow1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
							predhigh1=as.double(matrix(0,nrow=npred0,ncol=ntimeAll)),
							PACKAGE = "frailtypack") #19 arguments
							
						predMatLow <- matrix(ans$predlow1,nrow=npred0,ncol=ntimeAll)
						predMatHigh <- matrix(ans$predhigh1,nrow=npred0,ncol=ntimeAll)
						
					#=================================================#
					#************Prediction pour evenement de deces****
					#=================================================#
					}else{
						######### Distribution gamma #########
						if(fit$logNormal==0){
							for (k in 1:nrow(data)){
								realisations <- NULL
								for (i in 1:MC.sample){
									vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
									vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
									pred <- 1-((1+theta.mc[i]*(-log(vect.survival.X)))/(1+theta.mc[i]*(-log(vect.survival.X.horizon))))**(1/theta.mc[i])
									realisations <- cbind(realisations,pred)
								}
								predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
								predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
							}
						
						######### AK: Distribution LogNormale #########
						}else{
							mat.survival.X.mc <- NULL
							mat.survival.X.horizon.mc <- NULL

							for(i in 1:MC.sample){
								mat.survival.X.samp <- NULL
								mat.survival.X.horizon.samp <- NULL
								
								for(k in 1:nrow(data)){
									# vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
									vect.survival.X.samp <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
									# vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
									vect.survival.X.horizon.samp <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])
									mat.survival.X.samp <- rbind(mat.survival.X.samp,vect.survival.X.samp)
									mat.survival.X.horizon.samp <- rbind(mat.survival.X.horizon.samp,vect.survival.X.horizon.samp)
								}
								
								mat.survival.X.mc <- rbind(mat.survival.X.mc,mat.survival.X.samp)
								mat.survival.X.horizon.mc <- rbind(mat.survival.X.horizon.mc,mat.survival.X.horizon.samp)
							}		
							
							ans <- .Fortran("predict_LogN_sha",
								as.integer(nrow(data)),
								as.double(mat.survival.X),
								as.double(mat.survival.X.horizon),
								as.double(expBX),
								as.double(fit$sigma2),
								pred=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
								as.integer(1),
								as.integer(ntimeAll),
								as.integer(MC.sample),
								as.double(sigma2.mc),
								as.double(mat.survival.X.mc),
								as.double(mat.survival.X.horizon.mc),
								as.double(expBX.mc),
								predlow1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
								predhigh1=as.double(matrix(0,nrow=nrow(data),ncol=ntimeAll)),
								PACKAGE = "frailtypack")#15 arguments

							predMatLow <- matrix(ans$predlow1,nrow=nrow(data),ncol=ntimeAll)
							predMatHigh <- matrix(ans$predhigh1,nrow=nrow(data),ncol=ntimeAll)
						}
					}
					
				###############################
				###Prediction conditionnelle###
				###############################
				}else{					
					# cluster <- as.integer(as.vector(cluster))
					for (k in 1:nrow(data)){								
						realisations <- NULL
						frailty.mc <- NULL						
						mi <- fit$n.eventsbygrp[cluster[k]]	
						
						res1 <- sum((-log(sapply(tt1[which(clusterfit==cluster[k])],survival,ObjFrailty=fit))) %*% expBX3[which(clusterfit==cluster[k])])
						for (i in 1:MC.sample){
							vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
							######### Distribution gamma #########
							if(fit$logNormal==0){
								# if (k == 1) 
								frailty.mc <- c(frailty.mc,rgamma(1,shape=mi+1/theta.mc[i],scale=1/(res1+1/theta.mc[i])))
								pred <- 1-(vect.survival.X.horizon/vect.survival.X)**frailty.mc[i]
							######### AK: Distribution LogNormale #########
							}else{
								# if (k==1){
									res<-.Fortran("frailpred_sha_nor_mc",
										as.integer(fit$npar),
										frail.out=as.double(0),
										as.double(sigma2.mc[i]),
										as.double(res1),
										as.integer(mi),
										PACKAGE = "frailtypack" ) #5 arguments
									frailty.mc[i] <- res$frail.out
								# }
								pred <- 1-(vect.survival.X.horizon/vect.survival.X)**exp(frailty.mc[i])
							}
							realisations <- cbind(realisations,pred)
						}
						predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
						predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
					}
				}
				
			#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
			#-*-*-*-*-Pour un modele de Cox-*-*-*-*-#
			#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*#
			}else{
				for (k in 1:nrow(data)){
					realisations <- NULL
					for (i in 1:MC.sample){
						vect.survival.X <- sapply(sequence,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
						vect.survival.X.horizon <- sapply(sequence2,FUN=survival.mc,ObjFrailty=fit,para1=para.mc[i,],para2=para.mc2[i,])**expBX.mc[k,i]
						pred <- 1-(vect.survival.X.horizon/vect.survival.X)
						realisations <- rbind(realisations,pred)
					}
					predMatLow <- rbind(predMatLow,apply(realisations,1,quantile,probs=0.025))
					predMatHigh <- rbind(predMatHigh,apply(realisations,1,quantile,probs=0.975))
				}
			}
		} # Fin du calcul des bornes de confiance	
		
		out <- NULL
		out$call <- match.call()
		out$name.fit <- match.call()[[2]]
		
		if (event == 'Recurrent'){
			if(conditional) out$npred <- length(uni.cluster)
			else out$npred <- npred0
			out$event <- "Recurrent"
			out$predtimerec <- predtimerectmp
		}else{
			out$npred <- nrow(data)	
			out$event <- "Terminal"
		}
		
		out$moving.window <- moving.window
		if (moving.window){
			out$x.time <- sequence2
			out$t <- predTime
		}else{
			out$x.time <- sequence
		}
		out$pred <- predMat
		# colnames(out$pred) <- c("times",rep(" ",dim(out$pred)[2]-1))			
		colnames(out$pred) <- paste("time=", out$x.time)
		
		if (out$event == "Terminal") rownames(out$pred) <- paste("ind",1:out$npred)
		else{
			if (conditional) rownames(out$pred) <- paste("ind",unique(cluster))
			else rownames(out$pred) <- paste(nameGrup,unique(cluster)[order(unique(cluster))])
		}
		
		out$icproba <- ICproba
		if (ICproba){
			out$predLow <- predMatLow
			out$predHigh <- predMatHigh
			# colnames(out$predLow) <- c("times",rep(" ",dim(out$predLow)[2]-1))
			colnames(out$predLow) <- paste("time=", out$x.time)
			rownames(out$predLow) <- paste("ind",1:out$npred)
			# colnames(out$predHigh) <- c("times",rep(" ",dim(out$predHigh)[2]-1))
			colnames(out$predHigh) <- paste("time=", out$x.time)
			rownames(out$predHigh) <- paste("ind",1:out$npred)
		}
		if (fit$Frailty) {
			if (conditional) out$type <- 'conditional'
			else out$type <- 'marginal'
		}
		out$window <- window
		if (conditional) out$group <- unique(cluster)
		cat("Predictions done for",out$npred,"subjects \n")
		class(out) <- "predFrailty"
	}
	out
} 