## load frailtypack
library(frailtypack)

## load data
data(readmission)

## Cox
Cox <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
	n.knots=10,kappa1=10000,kappa2=10000,data=readmission,Frailty=FALSE)

## output
Cox

## summary
summary(Cox)	

## predictions
Cox$martingale.res
Cox$linear.pred

## plot
plot(Cox)

## Shared
sha <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
	n.knots=10,kappa1=10000,kappa2=10000,data=readmission,Frailty=TRUE)
	
## output
sha

## summary
summary(sha)

## predictions
sha$frailty.pred
sha$frailty.var
sha$martingale.res
sha$linear.pred

## plot
plot(sha)