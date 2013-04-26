
Cmeasures <- function(fitc, ties=1, marginal=0, cindex=0, Nboot=0){

	if(!(ties %in% c(0,1)))stop("Argument ties must be binary variable")
	if(!(marginal %in% c(0,1)))stop("Argument marginal must be binary variable")	
	if(!(cindex %in% c(0,1)))stop("Argument cindex must be binary variable")
	if((class(Nboot)!="numeric")||(Nboot < 0)||(Nboot > 1000)) stop("Argument Nboot must be a positive integer up to 1000")
	if(Nboot == 1)stop("More than one iteration is needed for bootstrap.")


	if(class(fitc)!="frailtyPenal")stop("the argument fitc must be from class frailtyPenal ")
	data <- eval(fitc$Names.data)
	
	update(fitc$formula,"~1")
	if(length(all.names(update(fitc$formula,"~1")))==4){
		Names.time <- as.character(fitc$formula[[2]][[2]])
		Names.event <- as.character(fitc$formula[[2]][[3]])
		surv.time <- data[,Names.time]
	}else{
		Names.time <- as.character(fitc$formula[[2]][[2]])
		Names.time2 <- as.character(fitc$formula[[2]][[3]])
		Names.event <- as.character(fitc$formula[[2]][[4]])
		surv.time <- data[,Names.time]-data[,Names.time2]

	}

	surv.status <- data[,Names.event]
	LPcond <- fitc$linear.pred
	if (fitc$Frailty){
		groupe <- sapply(data[,fitc$Names.cluster], FUN=function(x) which(sort(unique(data[,fitc$Names.cluster]))==x))
		LPmarg <- LPcond - fitc$frailty.pred[groupe] 
	}else{
		groupe <- seq(1,length(surv.status),1)
		LPmarg <- LPcond
	}
		
	BAW   <- cindexes.frailty(LPcond,LPmarg,surv.time,surv.status,groupe,ties,marginal,cindex)
	
	if(Nboot > 0){
		bootres <- boot(data=unique(groupe),statistic=statFP,R=Nboot, 
		fit=fitc, dataset=data, groupe=groupe, stimeboot=surv.time, statusboot=surv.status, ties=ties, marginal=marginal, cindex=cindex) 

		bootresCI <-apply(bootres$t[bootres$t[,1]==1,-1], MARGIN=2,FUN=function(x) quantile(x, probs=c(0.025,0.975), na.rm=TRUE))
		bootresSE <-apply(bootres$t[bootres$t[,1]==1,-1], MARGIN=2,FUN=function(x) sqrt(var(x)))
	}
	
	out <- NULL
	out$call <- fitc$formula
	out$Frailty <- fitc$Frailty	
	out$frequencies <- data.frame("Number.patients"=length(LPcond),"Number.events"=sum(surv.status),"Number.groups"=length(fitc$frailty.pred))
	out$Nboot <- Nboot
	if(Nboot > 0) out$Nbproblem <- sum(bootres$t[,1]!=1)
	out$ties <- ties

	out$CPEcond <- matrix(c(BAW$CPE.B.C,BAW$Npairs.between,BAW$CPE.W.C,BAW$Npairs.within,BAW$CPE.O.C,BAW$Npairs),nrow=2)
	if(Nboot > 0) out$CPEcond <- rbind(out$CPEcond,bootresSE[1:3],bootresCI[,1:3])
	colnames(out$CPEcond) <- c("Between","Within","Overall")


	if((ties==0) & (Nboot == 0))Namespair <- c("Concordance","Number pairs (no ties)")
	if((ties==0) & (Nboot > 0))Namespair <- c("Concordance","Number pairs (no ties)","SE","IC low","IC high")
	if((ties==1) & (Nboot == 0))Namespair <- c("Concordance","Number pairs")
	if((ties==1) & (Nboot > 0))Namespair <- c("Concordance","Number pairs","SE","IC low","IC high")

	rownames(out$CPEcond) <- Namespair

	if (!fitc$Frailty) {out$CPEcond <- t(t(out$CPEcond[,3])); colnames(out$CPEcond) <- "Overall"}	

	out$marginal <- marginal
	if(marginal==1){
		out$CPEmarg <- matrix(c(BAW$CPE.B.M,BAW$Npairs.between,BAW$CPE.W.M,BAW$Npairs.within,BAW$CPE.O.M,BAW$Npairs),ncol=3)
		if(cindex==1){
			if(Nboot > 0) out$CPEmarg <- rbind(out$CPEmarg,bootresSE[7:9],bootresCI[,7:9])
			
		}else{
			if(Nboot > 0) out$CPEmarg <- rbind(out$CPEmarg,bootresSE[4:6],bootresCI[,4:6])
		}

		colnames(out$CPEmarg) <- c("Between","Within","Overall")
		rownames(out$CPEmarg) <- Namespair
		if (!fitc$Frailty) {out$CPEmarg <- t(t(out$CPEmarg[,3])); colnames(out$CPEmarg) <- "Overall"}
	}

	out$cindex <- cindex	
	if(cindex==1){
		out$cindexcond <- matrix(c(BAW$cindex.B.C,BAW$comparable.between,BAW$cindex.W.C,BAW$comparable.within,BAW$cindex.O.C,BAW$comparable),ncol=3)
		if(Nboot > 0) out$cindexcond <- rbind(out$cindexcond,bootresSE[4:6],bootresCI[,4:6])
		colnames(out$cindexcond) <- c("Between","Within","Overall")
		rownames(out$cindexcond) <- Namespair

		if (!fitc$Frailty) {out$cindexcond <- t(t(out$cindexcond[,3])); colnames(out$cindexcond) <- "Overall"}	
	}


	if((marginal==1) & (cindex==1)){
		out$cindexmarg <- matrix(c(BAW$cindex.B.M,BAW$comparable.between,BAW$cindex.W.M,BAW$comparable.within,BAW$cindex.O.M,BAW$comparable),ncol=3)
		if(Nboot > 0) out$cindexmarg <- rbind(out$cindexmarg,bootresSE[10:12],bootresCI[,10:12])
		colnames(out$cindexmarg) <- c("Between","Within","Overall")
		rownames(out$cindexmarg) <- Namespair
		if (!fitc$Frailty) {out$cindexmarg <- t(t(out$cindexmarg[,3])); colnames(out$cindexmarg) <- "Overall"}	
	}	
	class(out) <- c("Cmeasures")
	out
  }
    


