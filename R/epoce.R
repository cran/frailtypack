
epoce <- function(fit, pred.times, newdata = NULL){

	if (missing(fit)) stop("The argument fit must be specified")
	if (class(fit)!="jointPenal") stop("The argument fit must be a class 'jointPenal'")
	if (missing(pred.times)) stop("The argument pred.times must be specified")
	if (class(pred.times)!="numeric") stop("pred.times must contain numerical values")

	if(!missing(newdata) & (class(newdata)!="data.frame")) stop("The argument newdata must be a 'data.frame'")

	nt <- length(pred.times)
	vopt <- fit$varHtotal
	b <- fit$b
	np <- length(fit$b)
	typeof <- fit$typeof
	nva <- fit$nvar

	if (typeof == 0){
		nz <- fit$n.knots
		zi <- fit$zi
	}else{
		nz <- 0
		zi <- 0
	}

	if (typeof == 1){
		nbintervR <- fit$nbintervR
		nbintervDC <- fit$nbintervDC
		ttt <- fit$time
		tttdc <- fit$timedc
	}else{
		nbintervR <- 0
		nbintervDC <- 0
		ttt <- 0
		tttdc <- 0
	}

	# recuperation des profils d'individus
	m <- fit$call
	m0 <- match.call()

	if (!missing(newdata)){
		if (length(colnames(eval(m$data)))!=length(colnames(eval(m0$newdata)))) stop("Your new dataset must have the same number of columns than the dataset used in the 'fit'")
		if (any(colnames(eval(m$data))!=colnames(eval(m0$newdata)))) stop("Your new dataset must have the very same variables than the dataset used in the 'fit'")
	}

	if (is.null(m$recurrentAG)) recurrentAG <- FALSE
	else recurrentAG <- TRUE

	m$formula.terminalEvent <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa1 <- m$kappa2 <- m$maxit <- m$hazard <- m$nb.int1 <-m$nb.int2 <- m$RandDist <- m$betaorder <- m$betaknots <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$init.Theta <- m$init.Alpha <- m$Alpha <- m$... <- NULL
	
	m[[1]] <- as.name("model.frame")
	if (!missing(newdata)) m[[3]] <- as.name(m0$newdata) # nouveau dataset

	dataset <- eval(m, sys.parent())

	typeofY <- attr(model.extract(dataset, "response"),"type")
	Y <- model.extract(dataset, "response")
	
	if (typeofY=="right"){
		tt0 <- rep(0,nobs)
		tt1 <- Y[,1]
		c <- Y[,2]
	}else{
		tt0 <- Y[,1]
		tt1 <- Y[,2]
		c <- Y[,3]
	}
	tt0 <- as.numeric(tt0)
	tt1 <- as.numeric(tt1)
	c <- as.numeric(c)

	class(m$formula) <- "formula"
	special <- c("strata", "cluster", "subcluster", "terminal", "num.id", "timedep")

	Terms <- terms(m$formula, special)#, data = m$data)
	
	m$formula <- Terms

	dropx <- NULL

	tempc <- untangle.specials(Terms, "cluster", 1:10)
	cluster <- strata(dataset[, tempc$vars], shortlabel = TRUE)
# 	numbers <- table(cluster)[order(unique(cluster))]
# 	newCluster <- rep(1:nsujet,numbers)
	dropx <- c(dropx,tempc$terms)

	tempterm <- untangle.specials(Terms, "terminal", 1:10)
	terminal <- strata(dataset[, tempterm$vars], shortlabel = TRUE)
	terminal <- as.numeric(as.character(terminal))
	dropx <- c(dropx,tempterm$terms)

	if (!is.null(dropx)) newTerms <- Terms[-dropx]
	else newTerms <- Terms

	X <- model.matrix(newTerms, dataset)
	if (ncol(X) > 1) X <- X[, -1, drop = FALSE]
	nva1 <- ncol(X)

	if (!missing(newdata)){
		nobs <- nrow(newdata)
		nsujet <- length(unique(cluster))
	}else{
		nobs <- fit$n
		nsujet <- fit$groups
	}

	if (!recurrentAG){
		tt0dc <- rep(0,nsujet)
		tt1dc <- aggregate(tt1,by=list(cluster),FUN=sum)[,2]
	}else{
		tt0dc <- rep(0,nsujet)
		tt1dc <- aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
	}
	cdc <- aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]

	m2 <- fit$call

	m2$n.knots <- m2$recurrentAG <- m2$cross.validation <- m2$kappa1 <- m2$kappa2 <- m2$maxit <- m2$hazard <- m2$nb.int1 <-m2$nb.int2 <- m2$RandDist <- m2$betaorder <- m2$betaknots <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$init.Theta <- m2$init.Alpha <- m2$Alpha <- m2$... <- NULL

	m2$formula[[3]] <- m2$formula.terminalEvent[[2]]
	m2$formula.terminalEvent <- NULL
	m2[[1]] <- as.name("model.frame")
	if (!missing(newdata)) m2[[3]] <- as.name(m0$newdata) # nouveau dataset
	datasetdc <- eval(m2, sys.parent())
	class(m2$formula) <- "formula"
	special2 <- c("strata", "timedep")
	Terms2 <- terms(m2$formula, special2)#, data = m3$data)

	X2 <- model.matrix(Terms2, datasetdc)
	if (ncol(X2) > 1) X2 <- X2[, -1, drop = FALSE]
	nva2 <- ncol(X2)

	if (!is.null(ncol(X2))){
		Xdc <- aggregate(X2[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]
		if (ncol(X2)>1){
			for (i in 2:ncol(X2)){
				Xdc.i <- aggregate(X2[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
				Xdc <- cbind(Xdc,Xdc.i)
			}
		}
	}else{
		Xdc <- aggregate(X2,by=list(cluster), FUN=function(x) x[length(x)])[,2]
	}

	if (!missing(newdata) & length(fit$coef)!=(ncol(X)+ncol(X2))) stop("Different covariates in model and newdata. Verify your dataset, be careful to the factor variables.")

	cat("\n")
	cat("Calculating ... \n")
	
	ans <- .Fortran("cvpl",
			as.integer(nobs),
			as.integer(nsujet),
			as.integer(cluster),
			as.integer(c),
			as.integer(cdc),
			as.integer(nva1),
			as.integer(nva2),
			as.double(X),
			as.double(Xdc),
			as.integer(typeof),
			as.integer(nz),
			as.double(zi),
			as.double(ttt),
			as.double(tttdc),
			as.integer(nbintervR),
			as.integer(nbintervDC),
			as.integer(np),
			as.double(b),
			as.double(vopt),
			as.double(tt0),
			as.double(tt1),
			as.double(tt0dc),
			as.double(tt1dc),
			as.integer(nt),
			as.double(pred.times),
			rl_cond=as.double(rep(0,nt)),
			epoir=as.double(rep(0,nt)),
			contribt=as.double(rep(0,nt*nsujet)),
			atrisk=as.double(rep(0,nt)),
			PACKAGE="frailtypack")

	out <- NULL
	if (!missing(newdata)) out$data <- m0$newdata
	else out$data <- fit$data
	out$new.data <- !is.null(newdata)
	out$pred.times <- pred.times
	out$mpl <- ans$rl_cond
	if (missing(newdata)) out$cvpl <- ans$epoir
	out$IndivContrib <- matrix(ans$contribt,nrow=nsujet,ncol=nt)
	out$AtRisk <- ans$atrisk

	cat("Estimators of EPOCE computed for",length(pred.times),"times \n")

	class(out) <- c("epoce")
	out
}
