

statFP <- function(data, indices, fit, dataset, groupe, stimeboot, statusboot, ties, marginal, cindex) {
	print("Bootstrap running ...")
	dataset$ligne <- 1:nrow(dataset)      
	databoot <- NA      
	for(i in 1:length(indices)) {
		datatmp<-cbind(dataset[groupe==indices[i],],groupeboot=i)
		databoot <- rbind(databoot,datatmp)
	}
	databoot <- databoot[! is.na(databoot$groupeboot),]
	databoot[,fit$Names.cluster] <- databoot$groupeboot
	indices.unit <- databoot$ligne
	model <- fit$call 
	model$data<- databoot
	fit.boot <- eval(model)
	istop <- fit.boot$istop
	LPcond <- fit.boot$linear.pred
	if (! is.null(fit.boot$frailty.pred)) LPmarg <- LPcond - fit.boot$frailty.pred[databoot$groupeboot] 
	else                                  LPmarg <- LPcond
	BAW.boot <- cindexes.frailty(LPcond, LPmarg, stimeboot[indices.unit], statusboot[indices.unit], databoot$groupe, ties,marginal,cindex)
	CPE.B.Cboot <- BAW.boot$CPE.B.C
	CPE.W.Cboot <- BAW.boot$CPE.W.C
	CPE.O.Cboot <- BAW.boot$CPE.O.C	
	if (cindex==1){
		cindex.B.Cboot <- BAW.boot$cindex.B.C
		cindex.W.Cboot <- BAW.boot$cindex.W.C
		cindex.O.Cboot <- BAW.boot$cindex.O.C      
	}     
	if (marginal==1){
		CPE.B.Mboot <- BAW.boot$CPE.B.M
		CPE.W.Mboot <- BAW.boot$CPE.W.M
		CPE.O.Mboot <- BAW.boot$CPE.O.M    
		
		if (cindex==1){
			cindex.B.Mboot <- BAW.boot$cindex.B.M
			cindex.W.Mboot <- BAW.boot$cindex.W.M
			cindex.O.Mboot <- BAW.boot$cindex.O.M      
		}
	}	
	res <- cbind(istop, CPE.B.Cboot, CPE.W.Cboot, CPE.O.Cboot)
	if (cindex==1)    res <- cbind(res, cindex.B.Cboot, cindex.W.Cboot, cindex.O.Cboot)
	if (marginal==1)  res <- cbind(res, CPE.B.Mboot, CPE.W.Mboot, CPE.O.Mboot)
	if (marginal==1 & cindex==1) res <- cbind(res, cindex.B.Mboot, cindex.W.Mboot, cindex.O.Mboot)
	return(res)
 }
