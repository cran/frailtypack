#' Identify individuals in Joint model for clustered data
#' 
#' This is a special function used in addition to the \code{cluster()} function
#' in the context of survival joint models for clustered data. This function
#' identifies subject index. It is used on the right hand side of a
#' 'frailtyPenal' formula. Using \code{num.id()} in a formula implies that a
#' joint frailty model for clustered data is fitted (Rondeau et al. 2011).
#' 
#' 
#' @usage num.id(x)
#' @param x A character or numeric variable which is supposed to indicate the
#' variable identifying individuals
#' @seealso \code{\link{frailtyPenal}}
#' @references V. Rondeau, J.P. Pignon, S. Michiels (2011). A joint model for
#' the dependence between clustered times to tumour progression and deaths: A
#' meta-analysis of chemotherapy in head and neck cancer. \emph{Statistical
#' methods in medical research} \bold{897}, 1-19.
#' @return No return value
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' \donttest{
#' 
#' data(readmission)
#' #-- here is generated cluster (5 clusters)
#' readmission <- transform(readmission,group=id%%5+1)
#' 
#' #-- exclusion all recurrent events --#
#' #--  to obtain framework of semi-competing risks --#
#' readmission2 <- subset(readmission, (t.start == 0 & event == 1) | event == 0)
#' 
#' joi.clus.gap <- frailtyPenal(Surv(time,event)~cluster(group)+
#' num.id(id)+dukes+charlson+sex+chemo+terminal(death),
#' formula.terminalEvent=~dukes+charlson+sex+chemo,
#' data=readmission2,recurrentAG=FALSE, n.knots=8,
#' kappa=c(1.e+10,1.e+10) ,Alpha="None")
#' 
#' }
#' 
#' 
"num.id"<-function(x)
 {
  x
 }
