\name{Cmeasures}
\alias{Cmeasures}
\alias{CbootstrapFP}
\alias{cindexes.frailty}
\alias{cindexes.W}
\alias{cindexes.B}
\alias{cindexes}
\alias{statFP}

\title{
Concordance measures in shared frailty and Cox models
}
\description{
Compute concordance probability estimation for Cox models or shared frailty models in case of grouped data (Mauguen et al. 2012). Concordance is given at different levels of comparison, taking into account the cluster membership: between-groups, within-groups and an overall measure, being a weighted average of the previous two. Can also compute the c-index (Harrell et al. 1996) at these three levels. It is possible to exclude tied pairs from concordance estimation (otherwise, account for 1/2). 
}
\usage{
Cmeasures(fitc, ties = 1, marginal = 0, cindex = 0, Nboot = 0)
}

\arguments{
  \item{fitc}{
A frailtyPenal object , for a shared frailty model. If the fit is a Cox model (Frailty=FALSE), no clustering membership is taken into account and only marginal concordance probability estimation is provided. Only an overall measure is given, where all patients are compared two by two.
If a counting process formulation is used to performed the fit, with 't.start' and 't.stop', the gap-times (t.stop-t.start) are used in the concordance estimation.

}
  \item{ties}{
Indicates if the tied pairs on prediction value must be included (ties=1) or excluded (ties=0) from the concordance estimation. Default is ties=1. When included, tied pairs account for 1/2 in the concordance.
}
  \item{marginal}{
Indicates if the concordance based on marginal predictions must be given (marginal=1) in addition to conditional ones or not (marginal=0). Marginal predictions do not include the frailty estimation in the linear predictor computation: uses "`Beta'X"' instead of "Beta'X + log z_i". Default is marginal=0.
}
  \item{cindex}{
Indicates if the c-index (Harrell et al. 1996) must be computed (cindex=1) in addition to the concordance probability estimation or not (cindex=0). C-index is also given at the three comparison levels (between, within and overall). Default is cindex=0.
}
  \item{Nboot}{
Number of bootstrap resamplings to compute standard-error of the concordances measures, as well as a percentile 95\% confidence interval. Nboot=0 indicates no bootstrap procedure. Maximum admitted is 1000. Minimum admitted is 2. Default is 0. Resampling is done at the group level. If Cox model is used, resampling is done at individual level.
}
}

\value{
\item{call}{The shared frailty model evaluated.}
\item{Frailty}{Logical value. Was model with frailties fitted.}
\item{frequencies}{Numbers of patients, events and groups used to fit the model.}
\item{Nboot}{Number of bootstrap resamplings required.}
\item{Nbproblem}{Number of bootstrap resamplings in which the model did not converge. Those are excluded when computing the 
			standard-error and the percentile 95\% confidence interval.}
\item{ties}{A binary, indicating if the tied pairs on prediction were used to compute the concordance.}
\item{CPEcond}{Values of the Concordance Probability Estimations (conditional) and the number of pairs used. If Nboot>0, give SE, the 
			standard-error of the parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper bounds of the percentile 
			confidence interval evaluated by bootstrap (2.5\% and 97.5\% percentiles).}
\item{marginal}{A binary, indicating if the marginal values were computed.}
\item{CPEmarg}{Values of the Concordance Probability Estimations (marginal)  and the number of pairs used, if marginal=1. If Nboot>0, 
			give SE, the standard-error of the parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper bounds of the 
			percentile confidence interval evaluated by bootstrap (2.5\% and 97.5\% percentiles).}
\item{cindex}{A binary, indicating if the c-indexes were computed.}
\item{cindexcond}{Values of the C-index (conditional)  and the number of pairs used. If Nboot>0, give SE, the standard-error of the 
			parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper bounds of the percentile confidence interval 
			evaluated by bootstrap (2.5\% and 97.5\% percentiles).}
\item{cindexmarg}{Values of the C-index (marginal) and the number of pairs used, if marginal=1. If Nboot>0, give SE, the 
			standard-error of the parameters evaluated by bootstrap, IC.low and IC.high, the lower and upper bounds of the percentile 
			confidence interval evaluated by bootstrap (2.5\% and 97.5\% percentiles).}
}
\references{
Concordance measures in shared frailty models: application to multicenter data in cancer prognosis. Mauguen A, Collette S, Pignon JP, Collette L and Rondeau V. (submitted 2012).

Harrell FE et al. Tutorial in biostatistics: multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Statistics in Medicine 1996; 15:361-387.

Gonen M, Heller G. Concordance probability and discriminatory power in proportional hazards regression. Biometrika 2005; 92:965-970.
}

\seealso{
\code{\link{print.Cmeasures}},
\code{\link{frailtyPenal}}
}
\examples{
## load data
data(readmission)

## a frailtypenal fit
fit.frailtypenal <-frailtyPenal(Surv(time,event)~as.factor(dukes)+cluster(id)+
as.factor(charlson)+as.factor(chemo),data=readmission,Frailty=TRUE,
cross.validation=FALSE,n.knots=10,kappa1=1,hazard="Splines")

## a Cmeasures call
fit.Cmeasures <- Cmeasures(fit.frailtypenal)
fit.Cmeasures.noties <- Cmeasures(fit.frailtypenal, ties=0)

## a short summary
fit.Cmeasures
fit.Cmeasures.noties
}

\keyword{concordance}