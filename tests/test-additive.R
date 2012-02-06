## load frailtypack
library(frailtypack)

## load data

data(dataAdditive) 

## fit 1

modAdd1 <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
correlation=FALSE,data=dataAdditive,n.knots=8,kappa1=10000)
		
## output
modAdd1

## summary

summary(modAdd1)

## predictions
modAdd1$frailty.pred
modAdd1$martingale.res
modAdd1$linear.pred

## plot
plot(modAdd1)


## fit 2

modAdd2 <- additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=10000)
		
## output
modAdd2

## summary

summary(modAdd2)

## predictions
modAdd2$frailty.pred
modAdd2$martingale.res
modAdd2$linear.pred

## plot
plot(modAdd2)



