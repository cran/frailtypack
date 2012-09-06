\name{hazard}
\alias{hazard}

\title{
Hazard function.
}
\description{
Let t be a continuous variable, we determine the value of the hazard function to t after run fit.
}
\usage{
hazard(t, ObjFrailty)
}

\arguments{
  \item{t}{
time for hazard function.
}
  \item{ObjFrailty}{
an object from the frailtypack fit.
}
}

\value{
return the value of hazard function in t.
}

\author{
V. Rondeau
}

\examples{

###  a fit Shared
data(readmission)
fit.shared <- frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+strata(sex),
             n.knots=10,kappa1=10000,kappa2=10000,data=readmission,
	     Frailty=TRUE)

### calling survival
hazard(20,fit.shared) 
}