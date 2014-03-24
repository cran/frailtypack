# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for nested model ...\n")

########################################################################
### NESTED frailty model
########################################################################

data("dataNested")

modNested <- frailtyPenal(Surv(t1, t2, event) ~ cluster(group) +
  subcluster(subgroup) + cov1 + cov2, data = dataNested,
  n.knots = 8, kappa1 = 50000, cross.validation = TRUE)

print(modNested, digits = 4)