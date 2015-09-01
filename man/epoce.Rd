
\name{epoce}
\alias{epoce}

\title{ Estimators of the Expected Prognostic Observed Cross-Entropy (EPOCE) for evaluating predictive accuracy of joint frailty models. }

\description{

This function computes estimators of the Expected Prognostic Observed Cross-Entropy (EPOCE) for evaluating the predictive accuracy of joint frailty models using \code{frailtyPenal}. On the same data as used for estimation of the joint frailty model, this function computes both the Mean Prognosis Observed Log-Likelihood (MPOL) and the Cross-Validated Prognosis Observed Log-Likelihood (CVPOL), two estimators of EPOCE. The latter corrects the MPOL estimate for over-optimism by approximated cross-validation. On external, this function only computes the Mean Prognostic Observed Log-Likelihood (MPOL).

}

\usage{
epoce(fit, pred.times, newdata)
}

\arguments{
  \item{fit}{A jointPenal object.}
  \item{pred.times}{Time or vector of times to compute epoce.}
  \item{newdata}{optional. When missing, the data used for estimating the fit are used, and CVPOL and MPOL are computed (internal validation). When newdata is specified, only MPOL is computed on this new dataset (external validation. The two datasets must have the same covariates with the same coding without missing data.}
}

\value{
  \item{data}{name of the data used to compute epoce}
  \item{new.data}{a boolean which is FALSE if computation is done on the same data as for estimation, and TRUE otherwise}
  \item{pred.times}{time or vector of times used in the function}
  \item{mpl}{values of MPOL for each pred.times}
  \item{cvpl}{values of CVPOL for each pred.times}
  \item{IndivContrib}{all the contributions to the log-likelihood for each pred.times}
  \item{AtRisk}{number of subject still at risk for each pred.times}
}

\references{
D. Commenges, B. Liquet, C. Proust-Lima (2012). Choice of prognostic estimators in joint models by estimating differences of expected conditional Kullback-Leibler risks. \emph{Biometrics}, \bold{68(2)}, 380-387.
}

\examples{

\dontrun{

########################################
#### epoce on a joint frailty model ####
########################################

data("readmission")

modJoint.gap <- frailtyPenal(Surv(t.start,t.stop,event)~ cluster(id) +
  dukes + charlson + sex + chemo + terminal(death),
  formula.terminalEvent = ~ dukes + charlson + sex + chemo ,
  data = readmission, n.knots = 8, kappa =c(2.11e+08,9.53e+11),
  recurrentAG=TRUE)

# computation on the same dataset
temps <- c(200,500,800,1100)
epoce <- epoce(modJoint.gap,temps)

print(epoce)
plot(epoce)

# computation on a new dataset, here a sample of readmission with the first 50 subjects
s <- readmission[1:100,]
epoce <- epoce(modJoint.gap,temps,newdata=s)

print(epoce)
plot(epoce)

}

}

\keyword{misc}






