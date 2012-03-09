## load frailtypack
library(frailtypack)

## load data
data(readmission)

## Shared
sha <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
	n.knots=10,kappa1=10000,kappa2=10000,data=readmission,Frailty=TRUE)
	
## Cmeasures
conc <- Cmeasures(sha,cindex=1,marginal=1)

## print
conc
