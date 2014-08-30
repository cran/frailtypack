
\name{Diffepoce}
\alias{Diffepoce}

\title{ Difference of Expected Prognostic Observed Cross-Entropy (EPOCE) estimators and its 95\% tracking interval between two joint frailty models. }

\description{

This function computes the difference of two EPOCE estimates (CVPOL and MPOL) and its 95\% tracking interval betwwen two joint frailty models estimated using \code{frailtyPenal}. Difference in CVPOL is computed when the EPOCE was previously estimated on the same dataset as used for estimation (using an approximated cross-validation), and difference in MPOL is computed when the EPOCE was previously estimated on an external datset.

}

\usage{
Diffepoce(epoce1, epoce2)
}

\arguments{
  \item{epoce1}{a first object inheriting from class epoce.}
  \item{epoce2}{a second object inheriting from class epoce.}
}

\value{
  \item{new.data}{a boolean which is FALSE if computation is done on the same data as for estimation, and TRUE otherwise}
  \item{pred.times}{time or vector of times used in the function}
  \item{DEPOCE}{the difference between the two MPOL or CVPOL for each time}
  \item{TIinf}{lower confidence band for the difference}
  \item{TIsup}{upper confidence band for the difference}
}

\details{

From the EPOCE estimates and the individual contributions to the prognostic observed log-likelihood obtained with \code{epoce} function on the same dataset from two different estimated joint frailty models, the difference of CVPOL (or MPOL) and its 95\% tracking interval is computed. The 95\% tracking interval is :
Delta(MPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for an external dataset
Delta(CVPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for the dataset used in \code{frailtyPenal}
where Delta(CVPOL) (or Delta(MPOL)) is the difference of CVPOL (or MPOL) of the two joint frailty models, and VARIANCE is the empirical variance of the difference of individuals contributions to the prognostic observed log-likelihoods of the two joint frailty models.

The estimators of EPOCE from  arguments epoce1 and epoce2 must have been computed on the same dataset and with the pred.times.

}

\references{
D. Commenges, B. Liquet, C. Proust-Lima (2012). Choice of prognostic estimators in joint models by estimating differences of expected conditional Kullback-Leibler risks. \emph{Biometrics}, \bold{68(2)}, 380-387.
}

\examples{

\dontrun{


data("readmission")

# first joint frailty model
joint1 <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
  dukes + charlson + sex + chemo + terminal(death),
  formula.terminalEvent = ~ dukes + charlson + sex + chemo ,
  data = readmission, n.knots = 8, kappa = c(2.11e+08,9.53e+11),
  recurrentAG=TRUE)

# second joint frailty model without dukes nor charlson as covariates
joint2 <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
  sex + chemo + terminal(death),
  formula.terminalEvent = ~ sex + chemo ,
  data = readmission, n.knots = 8, kappa = c(2.11e+08,9.53e+11),
  recurrentAG=TRUE)

temps <- c(200,500,800,1100)

# computation of estimators of EPOCE for the two models
epoce1 <- epoce(joint1,temps)
epoce2 <- epoce(joint2,temps)

# computation of the difference
diff <- Diffepoce(epoce1,epoce2)

print(diff)
plot(diff)

}

}

\keyword{misc}