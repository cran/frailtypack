#
#  names without 'factor' or 'as.factor'
#

factor.names <- function(x){
	x <- matrix(x,nrow=1)
	Names <- apply(x,MARGIN=2,FUN=function(x){
		if(length(grep("factor",x))!= 0){
			pos1 <- grep("\\(",unlist(strsplit(x,split="")))+1
			pos2 <- grep("\\)",unlist(strsplit(x,split="")))-1
			compris.factor <- substr(x,start=pos1,stop=pos2)
			after.factor <- substr(x,start=(pos2+2),stop=length(unlist(strsplit(x,split=""))))
			paste(compris.factor,after.factor,sep=".")
		}else{
			x
		}
	}
	)
	return(Names)
}




