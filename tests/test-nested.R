## load frailtypack
library(frailtypack)

## load data

data(dataNested)

## modClu

modClu<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,Frailty=TRUE,data=dataNested,
n.knots=8,kappa1=50000)

## output
modClu

## summary
summary(modClu)

## predictions
modClu$frailty.pred.group
modClu$frailty.pred.subgroup
modClu$martingale.res
modClu$linear.pred

## plot
plot(modClu)



