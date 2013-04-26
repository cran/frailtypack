
             
cindexes.frailty  <- function(lpcond,lpmarg,stime,status,groupe,ties,marginal,cindex) {

	cindexes.B.C  <- cindexes.B(lp=lpcond, stime=stime, status=status, groupe=groupe, ties=ties, cindex=cindex)
	cindexes.W.C  <- cindexes.W(lp=lpcond, stime=stime, status=status, groupe=groupe, ties=ties, cindex=cindex)
	
	if(marginal==1){
		cindexes.B.M  <- cindexes.B(lp=lpmarg, stime=stime, status=status, groupe=groupe, ties=ties, cindex=cindex)
		cindexes.W.M  <- cindexes.W(lp=lpmarg, stime=stime, status=status, groupe=groupe, ties=ties, cindex=cindex)
	}

	if(cindex==1){
		if (ties==1) {
			nwithincomp <-  cindexes.W.C$comparable
			nbetweencomp <- cindexes.B.C$comparable
		}
		if (ties==0) {
			nwithincomp <-  cindexes.W.C$comparable-cindexes.W.C$tiedcomp
			nbetweencomp <- cindexes.B.C$comparable-cindexes.B.C$tiedcomp
		}
		ncomp <- nwithincomp + nbetweencomp
	}
	
	if (ties==1) {    
		nwithintot <-  cindexes.W.C$Npairs
		nbetweentot <- cindexes.B.C$Npairs
	}
	if (ties==0) {    
		nwithintot <-  cindexes.W.C$Npairs-cindexes.W.C$tiedtot
		nbetweentot <- cindexes.B.C$Npairs-cindexes.B.C$tiedtot
	}
	ntot <- nwithintot + nbetweentot

	if(cindex==1) CVOAL.O.C  <- sum((nbetweencomp/ncomp)*cindexes.B.C$cindex, (nwithincomp/ncomp)*cindexes.W.C$cindex, na.rm=TRUE)
	CPE.O.C  <- sum((nbetweentot/ntot)*cindexes.B.C$CPE, (nwithintot/ntot)*cindexes.W.C$CPE, na.rm=TRUE)
	
	if(marginal==1){
		if(cindex==1) CVOAL.O.M  <- sum((nbetweencomp/ncomp)*cindexes.B.M$cindex, (nwithincomp/ncomp)*cindexes.W.M$cindex, na.rm=TRUE)
		CPE.O.M  <- sum((nbetweentot/ntot)*cindexes.B.M$CPE, (nwithintot/ntot)*cindexes.W.M$CPE, na.rm=TRUE)
	}

	out <- list(CPE.B.C=cindexes.B.C$CPE, CPE.W.C=cindexes.W.C$CPE, CPE.O.C=CPE.O.C,Npairs=ntot,Npairs.within=nwithintot, Npairs.between=nbetweentot)
	

	if(marginal==1){
		out <- c(out,CPE.B.M=cindexes.B.M$CPE, CPE.W.M=cindexes.W.M$CPE, CPE.O.M=CPE.O.M)
	}
	
	if(cindex==1){
		out <-c(out,cindex.B.C=cindexes.B.C$cindex, cindex.W.C=cindexes.W.C$cindex, cindex.O.C=CVOAL.O.C,comparable=ncomp, comparable.within=nwithincomp, comparable.between=nbetweencomp)
	}
	if((cindex==1) & (marginal==1)){
		out <-c(out,cindex.B.M=cindexes.B.M$cindex, cindex.W.M=cindexes.W.M$cindex, cindex.O.M=CVOAL.O.M)
	}
	return(out)
}
