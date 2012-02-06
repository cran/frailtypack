## load frailtypack
library(frailtypack)

## load data
data(readmission)

## modJoint_gap

modJoint_gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+as.factor(dukes)
	      +as.factor(charlson)+terminal(death),
	      formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
	      data=readmission,n.knots=14,kappa1=9550000000,kappa2=1410000000000,
	      Frailty=TRUE,joint=TRUE,recurrentAG=FALSE)
## output
modJoint_gap

## summary
summary(modJoint_gap)

## predictions
modJoint_gap$frailty.pred
modJoint_gap$frailty.var
modJoint_gap$martingale.res
modJoint_gap$martingaledeath.res
modJoint_gap$linear.pred
modJoint_gap$lineardeath.pred

## plot
plot(modJoint_gap)



