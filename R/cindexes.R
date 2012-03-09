

cindexes  <- function(lp,stime,status,ties,cindex) { 
	n <- length(lp);
	sumCPE <- 0;
	comparable  <- unusable  <- concordante  <- discordante  <- tiedcomp <- tiedtot  <- tiedtime  <- sumCPEj <- rep(NA,(n-1)) 
	Npairs <- n*(n-1)/2
	
	for (i in 1:(n-1)) {
		if(cindex==1){
			tiedtimej     <- (stime[i]==stime & status[i]==1 & status==1)
			tiedtime[i]   <- sum(tiedtimej[(i+1):n])
			unusablej     <- tiedtimej | ((stime[i]<stime & status[i]==0) | (stime<stime[i] & status==0)
						| (stime[i]==stime & status[i]==0 & status==0))
			unusable[i]   <- sum(unusablej[(i+1):n])
			comparablej   <- ! unusablej
			comparable[i] <- sum(comparablej[(i+1):n])
			concordantej  <- (comparablej & ((stime[i]<stime & lp[i]>lp) | (stime<stime[i] & lp>lp[i])
					| (stime[i]==stime & status[i]==0 & status==1 & lp>lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp<lp[i])))
			concordante[i]<- sum(concordantej[(i+1):n])
			discordantej  <- (comparablej & ((stime[i]>stime & lp[i]>lp) | (stime>stime[i] & lp>lp[i]) 
					| (stime[i]==stime & status[i]==0 & status==1 & lp<lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp>lp[i])))
			discordante[i]<- sum(discordantej[(i+1):n])
			tiedcompj         <- (comparablej & lp[i]==lp)
			tiedcomp[i]       <- sum(tiedcompj[(i+1):n])  
		}
		tiedtotj         <- (lp[i]==lp)
		tiedtot[i]       <- sum(tiedtotj[(i+1):n])
	
		bxjxi <- lp - lp[i]
		bxixj <- - bxjxi
		if (ties==1) sumCPEj[i] <- sum(((bxjxi<=0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj)))[(i+1):n])
		if (ties==0) sumCPEj[i] <- sum(((bxjxi<0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj)))[(i+1):n])
	}
	if (ties==1) {
		if(cindex==1) cindex_global <- (sum(concordante + tiedcomp/2))/sum(comparable)
		res.cpe <- (2/(n*(n-1))) * sum(sumCPEj)
	}
	if (ties==0) {
		if(cindex==1) cindex_global <- sum(concordante)/(sum(comparable)-sum(tiedcomp))
		res.cpe <- (1/((n*(n-1)/2)-sum(tiedtot))) * sum(sumCPEj)
	}

	out <- list(CPE=res.cpe,Npairs=Npairs,tiedtot=sum(tiedtot))
	if(cindex==1){
		out <- c(out,comparable=sum(comparable),concordante=sum(concordante),discordante=sum(discordante),
		tiedcomp=sum(tiedcomp),tiedtime=sum(tiedtime),unusable=sum(unusable),cindex=cindex_global)
	
	}
	return(out)
	
}
