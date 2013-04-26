

CbootstrapFP <- function(Nboot, dataset, fit, stimeboot, statusboot, groupe, frail, n.knots, kappa1) {

	bootres <- boot(data=unique(groupe),statistic=statFP,R=Nboot, 
	fit=fit, dataset=dataset, groupe=groupe, stimeboot=stimeboot, statusboot=statusboot, frail=frail,
	n.knots=n.knots, kappa1=kappa1) 

	bootresCI <-apply(bootres$t, MARGIN=2,FUN=function(x) quantile(x, probs=c(0.025,0.975)))
	bootresSE <-apply(bootres$t, MARGIN=2,FUN=function(x) sqrt(var(x)))
     
     return(list(bootres$t0, bootres$t, bootresCI, bootresSE))
}

