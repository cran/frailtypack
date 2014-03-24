# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for shared model ...\n")


########################################################################
### COX proportionnal hazard model with gap times
########################################################################

data("readmission")

mod.cox.gap <- frailtyPenal(Surv(time, event) ~ dukes + 
  charlson + sex + chemo, n.knots = 10, kappa1 = 1,
  data = readmission, cross.validation = TRUE)

print(mod.cox.gap, digits = 4)

########################################################################
### Shared frailty model with gap times
########################################################################

mod.sha.gap <- frailtyPenal(Surv(time,event) ~ cluster(id) +
  dukes + charlson + sex + chemo, n.knots = 10, kappa1 = 1,
  data = readmission, cross.validation = TRUE)
  
print(mod.sha.gap, digits = 4)

#########################################################################
### Stratified shared frailty model with gap times
#########################################################################

mod.sha.str.gap <- frailtyPenal(Surv(time, event) ~ cluster(id) +
  charlson + dukes + chemo + strata(sex), n.knots = 10, 
  kappa1 = 2.11e+08, kappa2 = 2.11e+08, data = readmission)

print(mod.sha.str.gap, digits = 4)


########################################################################
### Figures
########################################################################

pdf(file = "fig1.pdf", height = 3.6, width = 8.1)
par(mfrow = c(1, 3))
plot(mod.cox.gap, type.plot = "survival", main = "Cox model", conf.bands = TRUE)
plot(mod.sha.gap, type.plot = "survival", main = "Shared", conf.bands = TRUE)
plot(mod.sha.str.gap, type.plot = "survival", main = "Shared + Stratification",
  conf.bands = TRUE, pos.legend = "bottomleft", cex.legend = 1)