#' Plot Method for a joint surrogate mediation analysis model.
#'
#' Plots the estimated functions associated with the mediation analysis, i.e.
#' \eqn{g(s)}, \eqn{R(t)}
#' as well as the natural direct, indirect and total effects.
#' An option to plot the confidence bands of the function \eqn{g(s)} is available.
#' This option is also implemented for the confidence bands of the functions
#' \eqn{R(t)} and of the natural effects if these confidence bands are available.
#'
#' @aliases plot.jointSurroMed
#' @usage
#' \method{plot}{jointSurroMed}(x,plot.mediation="All",type.plot="Hazard",
#'	conf.bands=TRUE,endpoint=2,
#'	legend.pos = "topleft",...)
#' @param x An object of class \code{jointSurroMed} from a joint surrogate model
#'					with a mediation analysis
#' for longitudinal outcome and a terminal event, i.e., an
#' output from calling \code{jointSurroPenal} function with the option
#' 'mediation' set to TRUE.
#' @param plot.mediation A character string specifying the desired plot.
#'  Possible values are "All", "g","Rt" or "Effects". The default is
#' "All" which displays all three plots.
#' @param type.plot A character string specifying the type of curve
#' for the baseline hazards functions. Possible
#' value are "Hazard", or "Survival".
#' @param endpoint An integer specifying for which endpoint should
#' the baseline curves be plotted. Possible values are 0
#' for the surrogate endpoint only and 1 for the final endpoint or 2 for both.
#' Default is 2.
#' @param conf.bands Logical value. Determines whether confidence bands should be
#' plotted. The default is to do so if the confidence bands are available.
#' @param legend.pos The location of the legend can be specified by setting
#' this argument to a single keyword from the list '"bottomright"', '"bottom"',
#' '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#' '"center"'. The default is '"topleft"'
#' @param ... other unused arguments.
#' @return Print one or several plots for the mediation analysis
#' of a joint surrogate model
#' @seealso \code{\link{jointSurroPenal}}
#' @keywords file
#' @export

"plot.jointSurroMed" <- function (x, plot.mediation="All",type.plot="Hazard", conf.bands=TRUE,endpoint=2,
																	legend.pos = "topleft",...)
{
		#options check
	  plot.type <- charmatch(plot.mediation, c("All", "g","Rt","Effects"),nomatch = 0)
	  if (plot.type == 0) {
	    stop("'plot.mediation' must be one of the following: 'All', 'g','Rt' or 'Effects'")
	  }
		if(!(type.plot %in% c("Hazard","Survival"))){
			stop("'type.plot' must be one of the following: 'Hazard' or 'Survival'.")
		}
		if(!(endpoint %in% c(0,1,2))){
			stop("'endpoint' should be one of the following: 0,1 or 2.")
		}
		if(!(legend.pos %in% c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center" ))){
			stop("'legend.pos should be one of the following: 'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'")
		}
		if(type.plot=="Hazard"){
			if(endpoint==0){
				if(conf.bands){
					yupp = ifelse(max(x$lamS[,3])<0,0.8*max(x$lamS[,3]),1.2*max(x$lamS[,3]))
					ylow = ifelse(min(x$lamS[,2])<0,1.2*max(x$lamS[,2]),0.8*min(x$lamS[,2]))
					yy=x$lamS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the surrogate endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$lamS[,2],type='l',...)
					lines(xx,x$lamS[,3],type='l',...)
				}else{
					yupp = ifelse(max(x$lamS[,1])<0,0.8*max(x$lamS[,1]),1.2*max(x$lamS[,1]))
					ylow = ifelse(min(x$lamS[,1])<0,1.2*min(x$lamS[,1]),0.8*min(x$lamS[,1]))
					yy=x$lamS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the surrogate endpoint.",...)
				}
			}else if(endpoint==1){
				if(conf.bands){
					yupp = ifelse(max(x$lamT[,3])<0,0.8*max(x$lamT[,3]),1.2*max(x$lamT[,3]))
					ylow = ifelse(min(x$lamT[,2])<0,1.2*min(x$lamT[,2]),0.8*min(x$lamT[,2]))
					yy=x$lamT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the final endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$lamT[,2],type='l',lty=3,...)
					lines(xx,x$lamT[,3],type='l',lty=3,...)
				}else{
					yupp = ifelse(max(x$lamT[,1])<0,0.8*max(x$lamT[,1]),1.2*max(x$lamT[,1]))
					ylow = ifelse(min(x$lamT[,1])<0,1.2*min(x$lamT[,1]),0.8*min(x$lamT[,1]))
					yy=x$lamT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the final endpoint.",...)
				}
			}else{
				if(conf.bands){
					## S first ...
					yupp = ifelse(max(x$lamS[,3])<0,0.8*max(x$lamS[,3]),1.2*max(x$lamS[,3]))
					ylow = ifelse(min(x$lamS[,2])<0,1.2*min(x$lamS[,2]),0.8*min(x$lamS[,2]))
					yy=x$lamS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the surrogate endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$lamS[,2],type='l',lty=3,...)
					lines(xx,x$lamS[,3],type='l',lty=3,...)
					## Then T
					invisible(readline(prompt="Press [Enter] to continue:"))
					yupp = ifelse(max(x$lamT[,3])<0,0.8*max(x$lamT[,3]),1.2*max(x$lamT[,3]))
					ylow = ifelse(min(x$lamT[,2])<0,1.2*min(x$lamT[,2]),0.8*min(x$lamT[,2]))
					yy=x$lamT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the final endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$lamT[,2],type='l',lty=3,...)
					lines(xx,x$lamT[,3],type='l',lty=3,...)
				}else{
					## S first ...
					yupp = ifelse(max(x$lamS[,1])<0,0.8*max(x$lamS[,1]),1.2*max(x$lamS[,1]))
					ylow = ifelse(min(x$lamS[,1])<0,1.2*min(x$lamS[,1]),0.8*min(x$lamS[,1]))
					yy=x$lamS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the surrogate endpoint.",...)
					## Then T
					invisible(readline(prompt="Press [Enter] to continue:"))
					yupp = ifelse(max(x$lamT[,1])<0,0.8*max(x$lamT[,1]),1.2*max(x$lamT[,1]))
					ylow = ifelse(min(x$lamT[,1])<0,1.2*min(x$lamT[,1]),0.8*min(x$lamT[,1]))
					yy=x$lamT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Hazard',
							 main="Estimated baseline hazard function for the final endpoint.",...)
				}
			}
		}else{
			if(endpoint==0){
				if(conf.bands){
					yupp = ifelse(max(x$survS[,3])<0,0.8*max(x$survS[,3]),1.2*max(x$survS[,3]))
					ylow = ifelse(min(x$survS[,2])<0,1.2*min(x$survS[,2]),0.8*min(x$survS[,2]))
					yy=x$survS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the surrogate endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$survS[,2],type='l',lty=3,...)
					lines(xx,x$survS[,3],type='l',lty=3,...)
				}else{
					yupp = ifelse(max(x$survS[,1])<0,0.8*max(x$survS[,1]),1.2*max(x$survS[,1]))
					ylow = ifelse(min(x$survS[,1])<0,1.2*min(x$survS[,1]),0.8*min(x$survS[,1]))
					yy=x$survS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the surrogate endpoint.",...)
				}
			}else if(endpoint==1){
				if(conf.bands){
					yupp = ifelse(max(x$survT[,2])<0,0.8*max(x$survT[,2]),1.2*max(x$survT[,2]))
					ylow = ifelse(min(x$survT[,3])<0,1.2*min(x$survT[,3]),0.8*min(x$survT[,3]))
					yy=x$survT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the final endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$survT[,2],type='l',lty=3,...)
					lines(xx,x$survT[,3],type='l',lty=3,...)
				}else{
					yupp = ifelse(max(x$survT[,1])<0,0.8*max(x$survT[,1]),1.2*max(x$survT[,1]))
					ylow = ifelse(min(x$survT[,1])<0,1.2*min(x$survT[,1]),0.8*min(x$survT[,1]))
					yy=x$survT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the final endpoint.",...)
				}
			}else{
				if(conf.bands){
					## S first ...
					yupp = ifelse(max(x$survS[,3])<0,0.8*max(x$survS[,3]),1.2*max(x$survS[,3]))
					ylow = ifelse(min(x$survS[,2])<0,1.2*min(x$survS[,2]),0.8*min(x$survS[,2]))
					yy=x$survS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the surrogate endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$survS[,2],type='l',lty=3,...)
					lines(xx,x$survS[,3],type='l',lty=3,...)
					## Then T
					invisible(readline(prompt="Press [Enter] to continue:"))
					yupp = ifelse(max(x$survT[,2])<0,0.8*max(x$survT[,2]),1.2*max(x$survT[,2]))
					ylow = ifelse(min(x$survT[,3])<0,1.2*min(x$survT[,3]),0.8*min(x$survT[,3]))
					yy=x$survT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="HERE Estimated baseline survival function for the final endpoint\nand its 95% confidence bands.",...)
					lines(xx,x$survT[,2],type='l',lty=3,...)
					lines(xx,x$survT[,3],type='l',lty=3,...)
				}else{
					## S first ...
					yupp = ifelse(max(x$survS[,1])<0,0.8*max(x$survS[,1]),1.2*max(x$survS[,1]))
					ylow = ifelse(min(x$survS[,1])<0,1.2*min(x$survS[,1]),0.8*min(x$survS[,1]))
					yy=x$survS[,1]
					xx=x$xS
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the surrogate endpoint.",...)
					## Then T
					invisible(readline(prompt="Press [Enter] to continue:"))
					yupp = ifelse(max(x$survT[,1])<0,0.8*max(x$survT[,1]),1.2*max(x$survT[,1]))
					ylow = ifelse(min(x$survT[,1])<0,1.2*min(x$survT[,1]),0.8*min(x$survT[,1]))
					yy=x$survT[,1]
					xx=x$xT
					plot(xx,yy,ylim=c(ylow,yupp),type='l',
							col="black", xlab="Time",ylab='Survival',
							 main="Estimated baseline survival function for the final endpoint.",...)
				}
			}
		}
		if(plot.mediation=="g"){
			x<-x$mediation
			if(conf.bands==TRUE){
				data.rt<-x$data.rt
				rt.ci=x$Rt.ci
				nie.ci=x$NIE.ci
				nde.ci = x$NDE.ci
				te.ci =x$TE.ci
				data.g=x$data.g

				invisible(readline(prompt="Press [Enter] to continue:"))
				# plot gamma
				ymin.g<-ifelse(min(data.g$lower,na.rm=T)<0,
										 1.2*min(data.g$lower),
										 0.8*min(data.g$lower))
				ymax.g<-ifelse(max(data.g$upper,na.rm = T)<0,
										 0.8*max(data.g$upper,na.rm=T),
										 1.2*max(data.g$upper,na.rm=T))
				plot(x<-data.g$s,y<-data.g$g,type='l',col="black",
						 xlab="Time",ylab='g(S)',
						 main="Estimated function g(s) with its 95% confidence bands",
						 ylim=c(ymin.g,ymax.g),...)
				lines(x<-data.g$s,y<-data.g$upper,type='l',lty=3,...)
				lines(x<-data.g$s,y<-data.g$lower,type='l',lty=3,...)
			}else{
				#plot without confidence bands
				data.rt<-x$data.rt
				data.g<-x$data.g
				# plot gamma(s)
				invisible(readline(prompt="Press [Enter] to continue:"))
				plot(x<-data.g$s,y<-data.g$g,type='l',col="black",
						 xlab="Time",ylab='g(S)',main="Estimated function g(s)",...)
			}
		}
		if(plot.mediation=="Rt"){
			x<-x$mediation
			if(length(x)==9 & conf.bands){
				data.rt<-x$data.rt
				rt.ci=x$Rt.ci
				nie.ci=x$NIE.ci
				nde.ci = x$NDE.ci
				te.ci =x$TE.ci
				data.g=x$data.g
				#plot r(t)
				ymin<-ifelse(min(rt.ci$lower,na.rm=T)<0,
										 1.2*min(rt.ci$lower,na.rm=T),
										 0.8*min(rt.ci$yuppower,na.rm=T))
				ymax<-ifelse(max(rt.ci$upper,na.rm = T)<0,
										 0.8*max(rt.ci$upper,na.rm=T),
										 1.2*max(rt.ci$upper,na.rm=T))
				invisible(readline(prompt="Press [Enter] to continue:"))
				plot(x<-data.rt$Time,y<-data.rt$Rt,type='l',col="black",
						 xlab="Time",ylab='R(t)',main="Estimated R(t) with its 95% confidence band",
						 ylim=c(ymin,ymax),...)
				lines(x<-data.rt$Time,y<-rt.ci$upper,type='l',lty=2,col="black",...)
				lines(x<-data.rt$Time,y<-rt.ci$lower,type='l',lty=2,col="black",...)
			}else{
				#plot without confidence bands
				data.rt<-x$data.rt
				data.g<-x$data.g
				#plot r(t)
				invisible(readline(prompt="Press [Enter] to continue:"))
				plot(x<-data.rt$Time,y<-data.rt$Rt,type='l',col="black",
						 xlab="Time",ylab='R(t)',main="Estimated R(t)",...)
			}
		}
		if(plot.mediation=="Effects"){
			x<-x$mediation
			if(length(x)==9 & conf.bands){
				data.rt<-x$data.rt
				rt.ci=x$Rt.ci
				nie.ci=x$NIE.ci
				nde.ci = x$NDE.ci
				te.ci =x$TE.ci
				data.g=x$data.g
				#plot effects
				invisible(readline(prompt="Press [Enter] to continue:"))
				miny <- min(te.ci$lower,nie.ci$lower,nde.ci$lower,na.rm=T)
				maxy <-max(te.ci$upper,nie.ci$upper,nde.ci$upper,na.rm=T)
				ymin<-ifelse(miny<0,1.2*miny,0.8*miny)
				ymax<-ifelse(maxy<0,0.8*maxy,1.2*maxy)

				plot(x<-data.rt$Time,y<-data.rt$TE,type='l',col="black",
						 xlab="Time",ylab='Effects',main="Estimated natural effects with their 95% confidence bands",
						 ylim=c(ymin,ymax),...)
				lines(x<-data.rt$Time,y<-data.rt$NDE,type='l',col="green",...)
				lines(x<-data.rt$Time,y<-data.rt$NIE,type='l',col="red",...)

				lines(x<-data.rt$Time,y<-te.ci$lower,type='l',lty=2,col="black",...)
				lines(x<-data.rt$Time,y<-nde.ci$lower,type='l',lty=2,col="green",...)
				lines(x<-data.rt$Time,y<-nie.ci$lower,type='l',lty=2,col="red",...)

				lines(x<-data.rt$Time,y<-te.ci$upper,type='l',lty=2,col="black",...)
				lines(x<-data.rt$Time,y<-nde.ci$upper,type='l',lty=2,col="green",...)
				lines(x<-data.rt$Time,y<-nie.ci$upper,type='l',lty=2,col="red",...)

				legend(legend.pos,legend=c("Total effect","Direct effect","Indirect effect"),
							 col=c("black","green","red"),lty=1)
			}else{
				#plot without confidence bands
				data.rt<-x$data.rt
				data.g<-x$data.g
				#plot effect
				invisible(readline(prompt="Press [Enter] to continue:"))

				ymin<-ifelse(min(data.rt$TE,data.rt$NIE,data.rt$TE)<0,
										 1.2*min(data.rt$TE,data.rt$NIE,data.rt$TE),
										 0.8*min(data.rt$TE,data.rt$NIE,data.rt$TE))

				ymax<-ifelse(max(data.rt$TE,data.rt$NIE,data.rt$TE)<0,
										 0.8*max(data.rt$TE,data.rt$NIE,data.rt$TE),
										 1.2*max(data.rt$TE,data.rt$NIE,data.rt$TE))

				plot(x<-data.rt$Time,y<-data.rt$TE,type='l',col="black",
						 xlab="Time",ylab='Effects',main="Estimated natural effects",
						 ylim=c(ymin,ymax),...)
				lines(x<-data.rt$Time,y<-data.rt$NIE,type='l',col="green",...)
				lines(x<-data.rt$Time,y<-data.rt$NDE,type='l',col="red",...)
				legend(legend.pos,legend=c("Total effect","Direct effect","Indirect effect"),
							 col=c("black","green","red"),lty=1)

			}
		}
	  if(plot.mediation=="All"){
			x<-x$mediation
	    if(length(x)==9 & conf.bands){
	      data.rt<-x$data.rt
	      rt.ci=x$Rt.ci
	      nie.ci=x$NIE.ci
	      nde.ci = x$NDE.ci
	      te.ci =x$TE.ci
	      data.g=x$data.g

				invisible(readline(prompt="Press [Enter] to continue:"))
	      # plot gamma
				ymin.g<-ifelse(min(data.g$lower,na.rm=T)<0,
										 1.2*min(data.g$lower),
										 0.8*min(data.g$lower))
				ymax.g<-ifelse(max(data.g$upper,na.rm = T)<0,
										 0.8*max(data.g$upper,na.rm=T),
										 1.2*max(data.g$upper,na.rm=T))
	      plot(x<-data.g$s,y<-data.g$g,type='l',col="black",
	           xlab="Time",ylab='g(s)',
						 main="Estimated function g(s) with 95% confidence bands",
					 	 ylim=c(ymin.g,ymax.g),...)
				lines(x<-data.g$s,y<-data.g$upper,type='l',lty=3,...)
				lines(x<-data.g$s,y<-data.g$lower,type='l',lty=3,...)

	      #plot r(t)
	      ymin<-ifelse(min(rt.ci$lower,na.rm=T)<0,
	                   1.2*min(rt.ci$lower,na.rm=T),
	                   0.8*min(rt.ci$lower,na.rm=T))
	      ymax<-ifelse(max(rt.ci$upper,na.rm = T)<0,
	                   0.8*max(rt.ci$upper,na.rm=T),
	                   1.2*max(rt.ci$upper,na.rm=T))
	      invisible(readline(prompt="Press [Enter] to continue:"))
	      plot(x<-data.rt$Time,y<-data.rt$Rt,type='l',col="black",
	           xlab="Time",ylab='R(t)',main="Estimated R(t) with its 95% confidence band",
	           ylim=c(ymin,ymax),...)
	      lines(x<-data.rt$Time,y<-rt.ci$upper,type='l',lty=2,col="black",...)
	      lines(x<-data.rt$Time,y<-rt.ci$lower,type='l',lty=2,col="black",...)

	      #plot effects
	      invisible(readline(prompt="Press [Enter] to continue:"))
	      miny <- min(te.ci$lower,nie.ci$lower,nde.ci$lower,na.rm=T)
	      maxy <-max(te.ci$upper,nie.ci$upper,nde.ci$upper,na.rm=T)
	      ymin<-ifelse(miny<0,1.2*miny,0.8*miny)
	      ymax<-ifelse(maxy<0,0.8*maxy,1.2*maxy)

	      plot(x<-data.rt$Time,y<-data.rt$TE,type='l',col="black",
	           xlab="Time",ylab='Effects',main="Estimated natural effects with their 95% confidence bands",
	           ylim=c(ymin,ymax),...)
	      lines(x<-data.rt$Time,y<-data.rt$NDE,type='l',col="green",...)
	      lines(x<-data.rt$Time,y<-data.rt$NIE,type='l',col="red",...)

	      lines(x<-data.rt$Time,y<-te.ci$lower,type='l',lty=2,col="black",...)
	      lines(x<-data.rt$Time,y<-nde.ci$lower,type='l',lty=2,col="green",...)
	      lines(x<-data.rt$Time,y<-nie.ci$lower,type='l',lty=2,col="red",...)

	      lines(x<-data.rt$Time,y<-te.ci$upper,type='l',lty=2,col="black",...)
	      lines(x<-data.rt$Time,y<-nde.ci$upper,type='l',lty=2,col="green",...)
	      lines(x<-data.rt$Time,y<-nie.ci$upper,type='l',lty=2,col="red",...)

	      legend(legend.pos,legend=c("Total effect","Direct effect","Indirect effect"),
	             col=c("black","green","red"),lty=1)
	    }else if(length(x)==5 & conf.bands){
			 data.rt<-x$data.rt
			 data.g=x$data.g

			 invisible(readline(prompt="Press [Enter] to continue:"))
			 # plot gamma
			 ymin.g<-ifelse(min(data.g$lower,na.rm=T)<0,
										1.2*min(data.g$lower),
										0.8*min(data.g$lower))
			 ymax.g<-ifelse(max(data.g$upper,na.rm = T)<0,
										0.8*max(data.g$upper,na.rm=T),
										1.2*max(data.g$upper,na.rm=T))
			 plot(x<-data.g$s,y<-data.g$g,type='l',col="black",
						xlab="Time",ylab='g(s)',
						main="Estimated function g(s) with its 95% confidence bands",
						ylim=c(ymin.g,ymax.g),...)
			 lines(x<-data.g$s,y<-data.g$upper,type='l',lty=3,...)
			 lines(x<-data.g$s,y<-data.g$lower,type='l',lty=3,...)

			 invisible(readline(prompt="Press [Enter] to continue:"))
			 plot(x<-data.rt$Time,y<-data.rt$Rt,type='l',col="black",
						xlab="Time",ylab='R(t)',main="Estimated R(t)",...)


			 #plot effect
			 invisible(readline(prompt="Press [Enter] to continue:"))

			 ymin<-ifelse(min(data.rt$TE,data.rt$NIE,data.rt$NDE)<0,
										1.2*min(data.rt$TE,data.rt$NIE,data.rt$NDE),
										0.8*min(data.rt$TE,data.rt$NIE,data.rt$NDE))

			 ymax<-ifelse(max(data.rt$TE,data.rt$NIE,data.rt$NDE)<0,
										0.8*max(data.rt$TE,data.rt$NIE,data.rt$NDE),
										1.2*max(data.rt$TE,data.rt$NIE,data.rt$NDE))

			 plot(x<-data.rt$Time,y<-data.rt$TE,type='l',col="black",
						xlab="Time",ylab='Effects',main="Estimated natural effects",
						ylim=c(ymin,ymax),...)
			 lines(x<-data.rt$Time,y<-data.rt$NIE,type='l',col="green",...)
			 lines(x<-data.rt$Time,y<-data.rt$NDE,type='l',col="red",...)
			 legend(legend.pos,legend=c("Total effect","Direct effect","Indirect effect"),
							col=c("black","green","red"),lty=1)


			}else{
	      #plot without confidence bands
	      data.rt<-x$data.rt
	      data.g<-x$data.g
	      # plot gamma(s)
				ymin.g<-ifelse(min(data.g$lower,na.rm=T)<0,
										 1.2*min(data.g$lower),
										 0.8*min(data.g$lower))
				ymax.g<-ifelse(max(data.g$upper,na.rm = T)<0,
										 0.8*max(data.g$upper,na.rm=T),
										 1.2*max(data.g$upper,na.rm=T))
				invisible(readline(prompt="Press [Enter] to continue:"))
	      plot(x<-data.g$s,y<-data.g$g,type='l',col="black",
	           xlab="Time",ylab='g(s)',
						 main="Estimated function g(s)",
					 	 ylim=c(ymin.g,ymax.g),...)
	      #plot r(t)
	      invisible(readline(prompt="Press [Enter] to continue:"))
	      plot(x<-data.rt$Time,y<-data.rt$Rt,type='l',col="black",
	           xlab="Time",ylab='R(t)',main="Estimated R(t)",...)


	      #plot effect
	      invisible(readline(prompt="Press [Enter] to continue:"))

	      ymin<-ifelse(min(data.rt$TE,data.rt$NIE,data.rt$NDE)<0,
	                   1.2*min(data.rt$TE,data.rt$NIE,data.rt$NDE),
	                   0.8*min(data.rt$TE,data.rt$NIE,data.rt$NDE))

	      ymax<-ifelse(max(data.rt$TE,data.rt$NIE,data.rt$NDE)<0,
	                   0.8*max(data.rt$TE,data.rt$NIE,data.rt$NDE),
	                   1.2*max(data.rt$TE,data.rt$NIE,data.rt$NDE))

	      plot(x<-data.rt$Time,y<-data.rt$TE,type='l',col="black",
	           xlab="Time",ylab='Effects',main="Estimated natural effects",
	           ylim=c(ymin,ymax),...)
	      lines(x<-data.rt$Time,y<-data.rt$NIE,type='l',col="green",...)
	      lines(x<-data.rt$Time,y<-data.rt$NDE,type='l',col="red",...)
	      legend(legend.pos,legend=c("Total effect","Direct effect","Indirect effect"),
	             col=c("black","green","red"),lty=1)

	    }
	  }
	  return(invisible())
}
