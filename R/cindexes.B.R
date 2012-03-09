
              
cindexes.B  <- function(lp, stime, status, groupe, ties, cindex) {

	n <- length(lp)
	Npairs  <- comparable  <- unusable  <- concordante  <- discordante  <- tiedcomp <- tiedtot  <- tiedtime  <- rep(NA,(n-1))
	sumCPEj <- rep(NA,(n-1));sumCPE <- 0;
  
	for (i in 1:(n-1)) {
  
		Npairs[i] <- sum((groupe[i] != groupe)[(i+1):n])
		if(cindex==1){	
			tiedtimej     <- (groupe[i] != groupe) & stime[i]==stime & status[i]==1 & status==1
			tiedtime[i]   <- sum(tiedtimej[(i+1):n])
			unusablej     <- (groupe[i] != groupe) & tiedtimej | ((stime[i]<stime & status[i]==0) 
					| (stime<stime[i] & status==0) | (stime[i]==stime & status[i]==0 & status==0))
			unusable[i]   <- sum(unusablej[(i+1):n])
			comparablej   <- (groupe[i] != groupe) & (! unusablej)
			comparable[i] <- sum(comparablej[(i+1):n])
			concordantej  <- (comparablej & ((stime[i]<stime & lp[i]>lp) | (stime<stime[i] & lp>lp[i])
					| (stime[i]==stime & status[i]==0 & status==1 & lp>lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp<lp[i])))
			concordante[i]<- sum(concordantej[(i+1):n])
			discordantej  <- (comparablej & ((stime[i]>stime & lp[i]>lp) | (stime>stime[i] & lp>lp[i]) 
					| (stime[i]==stime & status[i]==0 & status==1 & lp<lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp>lp[i])))
			discordante[i]<- sum(discordantej[(i+1):n])
			tiedcompj     <- (comparablej & lp[i]==lp)
			tiedcomp[i]   <- sum(tiedcompj[(i+1):n])  
		}  
			tiedtotj      <- (groupe[i] != groupe) & (lp[i]==lp)
			tiedtot[i]    <- sum(tiedtotj[(i+1):n])

		bxjxi <- lp - lp[i]
		bxixj <- - bxjxi
		if (ties==1) {sumCPEj <- ( ((bxjxi<=0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj))) * (groupe != groupe[i]) )}
		if (ties==0) {sumCPEj <- ( ((bxjxi<0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj))) * (groupe != groupe[i]) )}
		sumCPE[i] <- sum(sumCPEj[(i+1):n])
	}

	if (ties==1) {
		if(cindex==1) cindex_global <- (sum(concordante) + sum(tiedcomp)/2)/sum(comparable)
		res.cpe <- sum(sumCPE)/sum(Npairs)
	}
	if (ties==0) {
		if(cindex==1) cindex_global <- sum(concordante)/(sum(comparable)-sum(tiedcomp))
		res.cpe <- sum(sumCPE)/(sum(Npairs)-sum(tiedtot))
	}
	out <- list(CPE=res.cpe,Npairs=sum(Npairs),tiedtot=sum(tiedtot)) 
	if(cindex==1){	
		out <- c(out,comparable=sum(comparable), concordante=sum(concordante),discordante=sum(discordante),tiedcomp=sum(tiedcomp),tiedtime=sum(tiedtime),
		unusable=sum(unusable),cindex=cindex_global)
	}
	return(out)
}
