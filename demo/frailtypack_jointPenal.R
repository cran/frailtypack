# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for joint model ...\n")

#########################################################################
### JOINT frailty model with gap times
#########################################################################

data("readmission")

modJoint.gap <- frailtyPenal(Surv(time,event)~ cluster(id) +
  dukes + charlson + sex + chemo + terminal(death),
  formula.terminalEvent = ~ dukes + charlson + sex + chemo,
  data = readmission, n.knots = 10, kappa1 = 100, kappa2 = 100)

print(modJoint.gap, digits = 4)

### Print the hazard ratios
summary(modJoint.gap, level = 0.95)

########################################################################
### Figures
########################################################################
pdf(file="fig3.pdf", height = 3.6, width = 8.1)
par(mfrow = c(1, 3))
plot(modJoint.gap, type.plot = "survival", event = "recurrent", main = "Recurrent",
  conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1.2, ylim = c(0, 1.2))
plot(modJoint.gap, type.plot = "survival", event = "terminal", main = "Terminal",
  conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1.2, ylim = c(0, 1.2))
plot(modJoint.gap, type.plot = "survival", event = "both", main = "Both",
  conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1, ylim = c(0, 1.2))