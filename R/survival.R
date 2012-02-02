survival <- function(t,ObjFrailty){

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
		out <- .Fortran("survival",as.double(t),as.double(the),as.double(the1),as.integer(nz+2),
		as.double(zi),survival=as.double(c(0,0)),as.integer(nst),PACKAGE = "frailtypack")	
		if(class(ObjFrailty) == "jointPenal"){
			res <- c(res,out$survival)
		}else{
			if(ObjFrailty$n.strat == 2){
				res <- c(res,out$survival)
			}else{
				res <- c(res,out$survival[1])
			}
		}
		return(res)
	}

	if (ObjFrailty$typeof == 1){
		res <- NULL
		if(class(ObjFrailty) == "jointPenal"){
			m <- ObjFrailty$nbintervR + ObjFrailty$nbintervDC
			b <- ObjFrailty$b[1:m]
			time <- ObjFrailty$time
			timedc <- ObjFrailty$timedc
			if((ObjFrailty$x1 > t) || (max(ObjFrailty$xSu1) < t)) stop(" Time exceeds the range allowed ")
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$xSu2) < t)) stop(" Time exceeds the range allowed ")
			out <- .Fortran("survivalj_cpm",as.double(t),as.double(b),as.integer(ObjFrailty$nbintervR),
			as.integer(ObjFrailty$nbintervDC),as.double(time),as.double(timedc),
			survival=as.double(c(0,0)),PACKAGE = "frailtypack")
			res <- c(res,out$survival)
		}else{	
			m <- ObjFrailty$n.strat*ObjFrailty$nbintervR
			b <- ObjFrailty$b[1:m]
			time <- ObjFrailty$time
			if((ObjFrailty$x1 > t) || (max(ObjFrailty$xSu1) < t)) stop(" Time exceeds the range allowed ")
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$xSu2) < t)) stop(" Time exceeds the range allowed ")
			out <- .Fortran("survival_cpm",as.double(t),as.double(b),
			as.integer(ObjFrailty$n.strat),as.integer(ObjFrailty$nbintervR),
			as.double(time),survival=as.double(c(0,0)),PACKAGE = "frailtypack")
	
			if(ObjFrailty$n.strat == 2){
				res <- c(res,out$survival)
			}else{
				res <- c(res,out$survival[1])	
			}
		}
		return(res)
	}


	if (ObjFrailty$typeof == 2){
		if(!t)stop(" Use only for time greater than 0")
		res <- NULL
		sc1 <- ObjFrailty$shape.weib[1]
		sh1 <- ObjFrailty$scale.weib[1]	
		if((ObjFrailty$x1 > t) || (max(ObjFrailty$xSu1) < t)) stop(" Time exceeds the range allowed ")
		res <- c(res,exp(-(t/sc1)^sh1))
		if(class(ObjFrailty) == "jointPenal"){
			sc1 <- ObjFrailty$shape.weib[2]
			sh1 <- ObjFrailty$scale.weib[2]	
			if((ObjFrailty$x2 > t) || (max(ObjFrailty$xSu2) < t)) stop(" Time exceeds the range allowed ")
			res <- c(res,exp(-(t/sc1)^sh1))
		}else{
			if(ObjFrailty$n.strat == 2){
				sc1 <- ObjFrailty$shape.weib[2]
				sh1 <- ObjFrailty$scale.weib[2]	
				if((ObjFrailty$x2 > t) || (max(ObjFrailty$xSu2) < t)) stop(" Time exceeds the range allowed ")
				res <- c(res,exp(-(t/sc1)^sh1))
			}
		}
		return(res)
	}







}
