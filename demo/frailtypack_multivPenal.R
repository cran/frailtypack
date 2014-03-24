# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for multivariate model ...\n")

########################################################################
### MULTIVARIATE frailty model with gap times
########################################################################

data("dataMultiv")

modMultiv <- multivPenal(Surv(TIMEGAP,INDICREC)~ cluster(PATIENT) + 
  v1 + v2 + event2(INDICMETA) + terminal(INDICDEATH),
  formula.Event2 =~ v1 + v2 + v3, formula.terminalEvent =~ v1,
  data = dataMultiv, hazard = "Weibull")

print(modMultiv, digits = 4)