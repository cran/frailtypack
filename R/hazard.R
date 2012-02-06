"hazard" <- function(t,ObjFrailty){

	if (ObjFrailty$typeof == 0){

		nz <- ObjFrailty$n.knots
		the <- ObjFrailty$b[1:(nz+2)] * ObjFrailty$b[1:(nz+2)]
		zi <- ObjFrailty$zi

		res <- NULL
		if((ObjFrailty$x1 > t) || (max(ObjFrailty$x1) < t)) stop(" Time exceeds the range allowed ")
		if(class(ObjFrailty) == "jointPenal"){
			nb1 <- nz+3
			nb2 <- 2*(nz+2)
			b <- ObjFrailty$b[nb1:nb2]
			the1 <- b * b	
			nst <- 2
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
		}else{
			if(ObjFrailty$n.strat == 2){
				nb1 <- nz+3
				nb2 <- 2*(nz+2)
				b <- ObjFrailty$b[nb1:nb2]
				the1 <- b * b	
				nst <- 2
				if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
			}else{
				the1 <- rep(0,(nz+6))
				nst <- 1
			}
		}
		
		out <- .Fortran("risque",as.double(t),as.double(the),as.double(the1),as.integer(nz+2),
		as.double(zi),risque=as.double(c(0,0)),as.integer(nst),PACKAGE = "frailtypack")
		
		if(class(ObjFrailty) == "jointPenal"){
			res <- c(res,out$risque)
		}else{
			if(ObjFrailty$n.strat == 2){
				res <- c(res,out$risque)
			}else{
				res <- c(res,out$risque[1])
			}
		}
		return(res)	
	}

	if (ObjFrailty$typeof == 1){
		res <- NULL
		if((ObjFrailty$x1 > t) || (max(ObjFrailty$x1) < t)) stop(" Time exceeds the range allowed ")
		x1 <- matrix(ObjFrailty$x1,nrow=3,ncol=ObjFrailty$nbintervR)
		x1 <- rbind(x1,rep(t,ObjFrailty$nbintervR))
		ind <- apply(x1,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
		res <- c(res,ObjFrailty$lam[(which(ind==1)*3-1),1])

		if(class(ObjFrailty) == "jointPenal"){
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
			x2 <- matrix(ObjFrailty$x2,nrow=3,ncol=ObjFrailty$nbintervDC)
			x2 <- rbind(x2,rep(t,ObjFrailty$nbintervDC))
			ind <- apply(x2,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
			res <- c(res,ObjFrailty$lam2[(which(ind==1)*3-1),1])
		}else{
# shared additive nested
			if(ObjFrailty$n.strat == 2){
				if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
				x2 <- matrix(ObjFrailty$x2,nrow=3,ncol=ObjFrailty$nbintervR)
				x2 <- rbind(x2,rep(t,ObjFrailty$nbintervR))
				ind <- apply(x2,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
				res <- c(res,ObjFrailty$lam2[(which(ind==1)*3-1),1])
			}
		}
		return(res)
	}


	if (ObjFrailty$typeof == 2){
		if(!t)stop(" Use only for time greater than 0")
		res <- NULL
		if((ObjFrailty$x1 > t) || (max(ObjFrailty$x1) < t)) stop(" Time exceeds the range allowed ")
		sc1 <- ObjFrailty$shape.weib[1]
		sh1 <- ObjFrailty$scale.weib[1]		
		res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))

		if(class(ObjFrailty) == "jointPenal"){
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
			sc1 <- ObjFrailty$shape.weib[2]
			sh1 <- ObjFrailty$scale.weib[2]	
			res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))

		}else{
			if(ObjFrailty$n.strat == 2){
				if((ObjFrailty$x2 > t) || (max(ObjFrailty$x2) < t)) stop(" Time exceeds the range allowed ")
				sc1 <- ObjFrailty$shape.weib[2]
				sh1 <- ObjFrailty$scale.weib[2]	
				res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
			}
		}
		return(res)
	}	
}


