#' Concordance measures in shared frailty and Cox proportional hazard models
#' 
#' Compute concordance probability estimation for Cox proportional hazard or
#' shared frailty models in case of grouped data (Mauguen et al. 2012).
#' Concordance is given at different levels of comparison, taking into account
#' the cluster membership: between-groups, within-groups and an overall
#' measure, being a weighted average of the previous two. Can also compute the
#' c-index (Harrell et al. 1996) at these three levels. It is possible to
#' exclude tied pairs from concordance estimation (otherwise, account for 1/2).
#' 
#' 
#' @aliases Cmeasures CbootstrapFP cindexes.frailty cindexes.W cindexes.B
#' cindexes statFP
#' @usage
#' 
#' Cmeasures(fitc, ties = 1, marginal = 0, cindex = 0, Nboot = 0, tau = 0,
#' data.val)
#' @param fitc A frailtyPenal object, for a shared frailty model. If the fit is
#' a Cox model, no clustering membership is taken into account and only
#' marginal concordance probability estimation is provided. Only an overall
#' measure is given, where all patients are compared two by two. If a counting
#' process formulation is used to performed the fit, with 't.start' and
#' 't.stop', the gap-times (t.stop-t.start) are used in the concordance
#' estimation.
#' @param ties Indicates if the tied pairs on prediction value must be included
#' (ties=1) or excluded (ties=0) from the concordance estimation. Default is
#' ties=1. When included, tied pairs account for 1/2 in the concordance.
#' @param marginal Indicates if the concordance based on marginal predictions
#' must be given (marginal=1) in addition to conditional ones or not
#' (marginal=0). Marginal predictions do not include the frailty estimation in
#' the linear predictor computation: uses "`Beta'X"' instead of "Beta'X + log
#' z_i". Default is marginal=0.
#' @param cindex Indicates if the c-index (Harrell et al. 1996) must be
#' computed (cindex=1) in addition to the concordance probability estimation or
#' not (cindex=0). C-index is also given at the three comparison levels
#' (between, within and overall). Default is cindex=0.
#' @param Nboot Number of bootstrap resamplings to compute standard-error of
#' the concordances measures, as well as a percentile 95\% confidence interval.
#' Nboot=0 indicates no bootstrap procedure. Maximum admitted is 1000. Minimum
#' admitted is 2. Default is 0. Resampling is done at the group level. If Cox
#' model is used, resampling is done at individual level.
#' @param tau Time used to limit the interval on which the concordance is
#' estimated. Note that the survival function for the underlying censoring time
#' distribution needs to be positive at tau. If tau=0, the maximum of the
#' observed event times is used. Default is tau=0.
#' @param data.val A dataframe. It is possible to specify a different dataset
#' than the one used in the model input in the argument 'fitc'. This new
#' dataset will be a validation population and the function will compute new
#' concordance measures from the parameters estimated on the development
#' population. In this case for conditional measures, the frailties are a
#' posteriori predicted. The two datasets must have the same covariates with
#' the same coding without missing data.
#' @return \item{call}{The shared frailty model evaluated.}
#' \item{Frailty}{Logical value. Was model with frailties fitted.}
#' \item{frequencies}{Numbers of patients, events and groups used to fit the
#' model.} \item{Npairs}{Number of pairs of subjects, between-groups,
#' within-groups and over all the population. If cindex=1, number of comparable
#' (useable) pairs also available.} \item{Nboot}{Number of bootstrap
#' resamplings required.} \item{ties}{A binary, indicating if the tied pairs on
#' prediction were used to compute the concordance.} \item{CPEcond}{Values of
#' Gonen & Heller's measure (conditional). If Nboot>0, give SE, the
#' standard-error of the parameters evaluated by bootstrap, IC.low and IC.high,
#' the lower and upper bounds of the percentile confidence interval evaluated
#' by bootstrap (2.5\% and 97.5\% percentiles).} \item{Cunocond}{Values of
#' Uno's measure (conditional). If Nboot>0, give SE, the standard-error of the
#' parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper
#' bounds of the percentile confidence interval evaluated by bootstrap (2.5\%
#' and 97.5\% percentiles).} \item{marginal}{A binary, indicating if the
#' marginal values were computed.} \item{CPEmarg}{Values of Gonen & Heller's
#' measure (marginal), if marginal=1. If Nboot>0, give SE, the standard-error
#' of the parameters evaluated by bootstrap, IC.low and IC.high, the lower and
#' upper bounds of the percentile confidence interval evaluated by bootstrap
#' (2.5\% and 97.5\% percentiles).} \item{Cunomarg}{Values of Uno's measure
#' (marginal), if marginal=1. If Nboot>0, give SE, the standard-error of the
#' parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper
#' bounds of the percentile confidence interval evaluated by bootstrap (2.5\%
#' and 97.5\% percentiles).} \item{cindex}{A binary, indicating if the
#' c-indexes were computed.}
#' 
#' \item{cindexcond}{Values of the C-index of Harrell (conditional). If
#' Nboot>0, give SE, the standard-error of the parameters evaluated by
#' bootstrap, IC.low and IC.high, the lower and upper bounds of the percentile
#' confidence interval evaluated by bootstrap (2.5\% and 97.5\% percentiles).}
#' \item{cindexmarg}{Values of the C-index of Harrell (marginal), if
#' marginal=1. If Nboot>0, give SE, the standard-error of the parameters
#' evaluated by bootstrap, IC.low and IC.high, the lower and upper bounds of
#' the percentile confidence interval evaluated by bootstrap (2.5\% and 97.5\%
#' percentiles).}
#' @seealso \code{\link{print.Cmeasures}},\code{\link{frailtyPenal}}
#' @references
#' 
#' Mauguen, A., Collette, S., Pignon, J. P. and Rondeau, V. (2013). Concordance
#' measures in shared frailty models: application to clustered data in cancer
#' prognosis. \emph{Statistics in Medicine} \bold{32}, 27, 4803-4820
#' 
#' Harrell, F.E. et al. (1996). Tutorial in biostatistics: multivariable
#' prognostic models: issues in developing models, evaluating assumptions and
#' adequacy, and measuring and reducing errors. \emph{Statistics in Medicine}
#' \bold{15}, 361-387.
#' 
#' Gonen, M., Heller, G. (2005). Concordance probability and discriminatory
#' power in proportional hazards regression. \emph{Biometrika} \bold{92},
#' 965-970.
#' @keywords concordance
#' @export
#' @examples
#' 
#' 
#' \donttest{
#' 
#' #-- load data
#' data(readmission)
#' 
#' #-- a frailtypenal fit
#' fit <- frailtyPenal(Surv(time,event)~cluster(id)+dukes+
#' charlson+chemo,data=readmission,cross.validation=FALSE,
#' n.knots=10,kappa=1,hazard="Splines")
#' 
#' #-- a Cmeasures call
#' fit.Cmeasures <- Cmeasures(fit)
#' fit.Cmeasures.noties <- Cmeasures(fit, ties=0)
#' fit.Cmeasures.marginal <- Cmeasures(fit, marginal=1)
#' fit.Cmeasures.cindex <- Cmeasures(fit, cindex=1)
#' 
#' #-- a short summary
#' fit.Cmeasures
#' fit.Cmeasures.noties
#' fit.Cmeasures.marginal
#' fit.Cmeasures.cindex
#' 
#' }
#' 
#' 
Cmeasures <- function(fitc, ties=1, marginal=0, cindex=0, Nboot=0, tau=0, data.val){

	if (missing(fitc)) stop("Need a fit")
	if (!inherits(fitc, "frailtyPenal")) stop("The argument fitc must be a frailtyPenal object")
	if(!(ties %in% c(0,1)))stop("Argument ties must be binary variable")
	if(!(marginal %in% c(0,1)))stop("Argument marginal must be binary variable")	
	if(!(cindex %in% c(0,1)))stop("Argument cindex must be binary variable")
	if(!inherits(Nboot, "numeric")||(Nboot < 0)||(Nboot > 1000)) stop("Argument Nboot must be a positive integer up to 1000")
	if(Nboot == 1)stop("More than one iteration is needed for bootstrap.")

	if (!missing(data.val)) data <- data.val
	else data <- eval(fitc$Names.data)
	
	# enlever les NA du data s'il y en a
	drop <- rep(TRUE,nrow(data))
	if (any(is.na(data))){
		drop <- complete.cases(data)
		warning("Missing observations in data ( ",length(drop==FALSE)," rows removed )")
		data <- data[drop,]
	}
	
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
	
	if(tau==0) tau=max(surv.time[surv.status==1])
	
	if (!missing(data.val)){ # validation externe sur data.val
		m <- fitc$call
		m2 <- match.call()
		m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
		
		m[[1]] <- as.name("model.frame")
		m[[3]] <- as.name(m2$data.val)
		
		if (fitc$Frailty){
			m$formula <- unlist(strsplit(deparse(m$formula)," "))
			clus <- grep("cluster",m$formula)
			if (clus==length(m$formula)) m$formula <- as.formula(paste(m$formula[-c(clus,max(which(m$formula=="+")))],collapse=""))
			else m$formula <- as.formula(paste(m$formula[-clus],collapse=""))
		}
		m$formula[[2]] <- NULL # pas besoin du Surv dans formula
		
		dataset <- eval(m, sys.parent())
		dataset <- dataset[drop,] # enlever les NA si y en a
		
		class(m$formula) <- "formula"
		special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")
		
		Terms <- terms(m$formula, special, data = data)
		
		X <- model.matrix(Terms, dataset)
		if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
		
		if (length(fitc$coef)!=ncol(X)) stop("Different covariates in model and data.val. Verify your dataset, be careful to the factor variables")
		
		if (fitc$Frailty){
			# calcul de la fragilite pour le nouveau dataset
			tt1 <- surv.time
			tt1[tt1>max(fitc$x1)] <- max(fitc$x1)
			res1 <- as.vector(-log(sapply(tt1,survival,ObjFrailty=fitc)) * exp(fitc$coef %*% t(X)))
			Hazg <- aggregate(res1,by=list(data[,fitc$Names.cluster]),FUN=sum)[,2]
			
			mg <- as.vector(table(data[surv.status==1,][,fitc$Names.cluster]))
			uk <- rep((1/fitc$theta+mg)/(1/fitc$theta+Hazg),times=table(data[,fitc$Names.cluster]))
			
			LPcond <- fitc$coef %*% t(X) + log(uk)
		}else{
			LPcond <- fitc$coef %*% t(X)
		}
		
		if (fitc$Frailty){
			groupe <- sapply(data[,fitc$Names.cluster], FUN=function(x) which(sort(unique(data[,fitc$Names.cluster]))==x))
			LPmarg <- LPcond - log(uk)
		}else{
			groupe <- seq(1,length(surv.status),1)
			LPmarg <- LPcond
		}
	}else{
		LPcond <- fitc$linear.pred
		
		if (fitc$Frailty){
			groupe <- sapply(data[,fitc$Names.cluster], FUN=function(x) which(sort(unique(data[,fitc$Names.cluster]))==x))
			LPmarg <- LPcond - log(fitc$frailty.pred[groupe])
		}else{
			groupe <- seq(1,length(surv.status),1)
			LPmarg <- LPcond
		}
	}
	
	BAW   <- cindexes.frailty(LPcond,LPmarg,surv.time,surv.status,groupe,ties,tau,marginal,cindex)
	
	if(Nboot > 0){
		bootres <- boot(data=unique(groupe),statistic=statFP,R=Nboot, 
		fit=fitc, dataset=data, LPcond=LPcond, LPmarg=LPmarg, groupe=groupe, stimeboot=surv.time, statusboot=surv.status, ties=ties, tau=tau, marginal=marginal, cindex=cindex) 
		
# 		bootresCI <-apply(bootres$t[bootres$t[,1]==1,-1], MARGIN=2,FUN=function(x) quantile(x, probs=c(0.025,0.975), na.rm=TRUE))
# 		bootresSE <-apply(bootres$t[bootres$t[,1]==1,-1], MARGIN=2,FUN=function(x) sqrt(var(x)))
		bootresCI <-apply(bootres$t, MARGIN=2,FUN=function(x) quantile(x, probs=c(0.025,0.975), na.rm=TRUE))
		bootresSE <-apply(bootres$t, MARGIN=2,FUN=function(x) sqrt(var(x)))
	}
	
	out <- NULL
	out$call <- fitc$formula
	out$Frailty <- fitc$Frailty
	out$frequencies <- data.frame("Number.patients"=length(LPcond),"Number.events"=sum(surv.status),"Number.groups"=length(unique(groupe)))
	out$Npairs <- data.frame("Between"=c(BAW$Npairs.between,BAW$comparable.between),"Within"=c(BAW$Npairs.within,BAW$comparable.within),"Overall"=c(BAW$Npairs,BAW$comparable))
	if (cindex==1) rownames(out$Npairs) <- c("Nb pairs","Nb pairs useable") else rownames(out$Npairs) <- "Nb pairs"
	out$Nboot <- Nboot
	#if(Nboot > 0) out$Nbproblem <- sum(bootres$t[,1]!=1)
	out$ties <- ties
	if (!missing(data.val)) out$Names.data <- match.call()$data.val
	else out$Names.data <- fitc$Names.data
	out$marginal <- marginal
	out$cindex <- cindex
	
	if(marginal==0){
		out$CPEcond <- matrix(c(BAW$CPE.B.C,BAW$CPE.W.C,BAW$CPE.O.C),nrow=1)
		out$Cunocond <- matrix(c(BAW$Cuno.B.C,BAW$Cuno.W.C,BAW$Cuno.O.C),nrow=1)
		if(Nboot > 0){
			out$CPEcond <- rbind(out$CPEcond,bootresSE[1:3],bootresCI[,1:3])
			out$Cunocond <- rbind(out$Cunocond,bootresSE[4:6],bootresCI[,4:6])
		}
		if (!fitc$Frailty) {
			out$CPEcond <- t(t(out$CPEcond[,3]))
			out$Cunocond <- t(t(out$Cunocond[,3]))
		}
	}else{ # if (marginal==1)
		out$CPEmarg <- matrix(c(BAW$CPE.B.M,BAW$CPE.W.M,BAW$CPE.O.M),ncol=3)
		out$Cunomarg <- matrix(c(BAW$Cuno.B.M,BAW$Cuno.W.M,BAW$Cuno.O.M),ncol=3)
		if(Nboot > 0){
			out$CPEmarg <- rbind(out$CPEmarg,bootresSE[1:3],bootresCI[,1:3])
			out$Cunomarg <- rbind(out$Cunomarg,bootresSE[4:6],bootresCI[,4:6])
		}
		if (!fitc$Frailty) {
			out$CPEmarg <- t(t(out$CPEmarg[,3]))
			out$Cunomarg <- t(t(out$Cunomarg[,3]))
		}
	}
	
	if((marginal==0) & (cindex==1)){
		out$cindexcond <- matrix(c(BAW$cindex.B.C,BAW$cindex.W.C,BAW$cindex.O.C),ncol=3)
		if(Nboot > 0) out$cindexcond <- rbind(out$cindexcond,bootresSE[7:9],bootresCI[,7:9])
		if (!fitc$Frailty) out$cindexcond <- t(t(out$cindexcond[,3]))
	}
	
	if((marginal==1) & (cindex==1)){
		out$cindexmarg <- matrix(c(BAW$cindex.B.M,BAW$cindex.W.M,BAW$cindex.O.M),ncol=3)
		if(Nboot > 0) out$cindexmarg <- rbind(out$cindexmarg,bootresSE[7:9],bootresCI[,7:9])
		if (!fitc$Frailty) out$cindexmarg <- t(t(out$cindexmarg[,3]))
	}
	class(out) <- c("Cmeasures")
	out
  }
    


