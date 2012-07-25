
\name{frailtypack-package}
\Rdversion{2.2-24}
\alias{frailtypack-package}
\alias{frailtypack}
\docType{package}

\title{General Frailty models using a semi-parametrical penalized likelihood estimation or a parametrical estimation}

\description{
 Frailtypack now fits several classes of frailty models using a penalized likelihood estimation on the hazard function but also a parametrical estimation. 1) A shared gamma frailty model and Cox proportional hazard model. Clustered and recurrent survival times can be studied (the Andersen-Gill(1982) approach has been implemented for recurrent events). An automatic choice of the smoothing parameter is possible using an approximated cross-validation procedure. 2) Additive frailty models for proportional hazard models with two correlated random effects (intercept random effect with random slope). 3) Nested frailty models for hierarchically clustered data (with 2 levels of clustering) by including two iid gamma random effects. 4) Joint frailty models in the context of joint modelling of recurrent events with terminal event. Prediction values are available. Left truncated, right-censored data, interval-censored data (only for Cox proportional hazard and shared frailty model) and strata (max=2) are allowed. The package includes concordance measures for Cox proportional hazards models and for shared frailty models.
}
\details{
  \tabular{ll}{
     Package: \tab frailtypack\cr
     Type: \tab Package\cr
     Version: \tab 2.2-24\cr
     Date: \tab 2012-07-24\cr
     License: \tab GPL (>= 2.0)\cr
     LazyLoad: \tab yes\cr
  }  
 Frailtypack now fits several classes of frailty models using a penalized likelihood estimation on the hazard function but also a parametrical estimation. 1) A shared gamma frailty model and Cox proportional hazard model. Clustered and recurrent survival times can be studied (the Andersen-Gill(1982) approach has been implemented for recurrent events). An automatic choice of the smoothing parameter is possible using an approximated cross-validation procedure. 2) Additive frailty models for proportional hazard models with two correlated random effects (intercept random effect with random slope). 3) Nested frailty models for hierarchically clustered data (with 2 levels of clustering) by including two iid gamma random effects. 4) Joint frailty models in the context of joint modelling of recurrent events with terminal event. Prediction values are available. Left truncated, right-censored data, interval-censored data (only for Cox proportional hazard and shared frailty model) and strata (max=2) are allowed. The package includes concordance measures for Cox proportional hazards models and for shared frailty models.
}

\author{Virginie Rondeau, Juan R. Gonzalez, Yassin Mazroui, Audrey Mauguen, Amadou Diakite and Alexandre Laurent}

\references{
V. Rondeau, Y. Mazroui and J. R. Gonzalez (2012). Frailtypack: An R package for the analysis of correlated survival data with frailty models using penalized likelihood estimation or parametrical estimation.
\emph{Journal of Statistical Software} \bold{47}, 1-28.

V. Rondeau and J. R. Gonzalez. Frailtypack: A computer program for the analysis of correlated failure time data using penalized likelihood estimation. \emph{Computer Methods and Programs in Biomedicine, 80(2)}: \bold{154-164}, 2005.

V. Rondeau, S. Michiels, B. Liquet, and J. P. Pignon (2008). Investigating trial and treatment heterogeneity in an individual patient data meta-analysis of survival data by mean of the maximum penalized likelihood approach. \emph{Statistics in Medecine}, \bold{27}, 1894-1910.

V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P. Soubeyran (2007). Joint frailty models for recurring events and death using maximum penalized likelihood estimation:application on cancer events. \emph{Biostatistics}, \bold{8},4, 708-721.

V. Rondeau, D. Commenges, and P. Joly (2003). Maximum penalized likelihood estimation in a gamma-frailty model. \emph{Lifetime Data Analysis} \bold{9}, 139-153.

D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. \emph{SIAM Journal of Applied Mathematics}, 431-441. 

V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum penalized likelihood estimation. \emph{Statistics in Medecine}, \bold{25}, 4036-4052.

}

\seealso{
\code{\link{frailtyPenal}},
\code{\link{additivePenal}},
}


\examples{
\dontrun{
### Additive model with 1 covariate ###

 data(dataAdditive) 
 modAdd<-additivePenal(Surv(t1,t2,event)~cluster(group)+var1+slope(var1),
                correlation=TRUE,data=dataAdditive,n.knots=8,kappa1=10000,
		hazard="Splines")
		

### Joint model (recurrent and terminal events) with 2 covariates ###
### on a simulated dataset ###

 data(readmission)
 modJoint_gap<-frailtyPenal(Surv(time,event)~cluster(id)+sex+as.factor(dukes)
	      +as.factor(charlson)+terminal(death),
	      formula.terminalEvent=~sex+as.factor(dukes)+as.factor(charlson),
	      data=readmission,n.knots=14,kappa1=9550000000,kappa2=1410000000000,
	      Frailty=TRUE,joint=TRUE,recurrentAG=FALSE,hazard="Splines")

 ### Nested model (or hierarchical model) with 2 covariates ###
 
 data(dataNested)
 modClu<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+
 subcluster(subgroup)+cov1+cov2,Frailty=TRUE,data=dataNested,
 n.knots=8,kappa1=50000,hazard="Splines") 

 ###  Semi-parametrical Shared model   ###
 
 data(readmission)
 frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
          as.factor(charlson)+cluster(id),data=readmission, Frailty=TRUE,
          n.knots=6,kappa1=5000,recurrentAG=TRUE,cross.validation=TRUE,
	  hazard="Splines")
	  
###  Parametrical Shared model   ###

 data(readmission)
 frailtyPenal(Surv(t.start,t.stop,event)~as.factor(sex)+as.factor(dukes)+
          as.factor(charlson)+cluster(id),data=readmission, Frailty=TRUE,
	  hazard="Piecewise-per",nb.int1=6)
}
}

\keyword{package}
