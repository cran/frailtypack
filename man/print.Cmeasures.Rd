\name{print.Cmeasures}
\alias{print.Cmeasures}

\title{
Print a short summary of results of Cmeasure function.
}
\description{
Print a short summary of results of the concordance measure estimated by the Cmeasure function.
}
\usage{
\method{print}{Cmeasures}(x, digits = 3, ...)
}

\arguments{
  \item{x}{
a Cmeasures object.
}
  \item{digits}{
Number of digits to print
}
  \item{\dots}{
Other unused arguments
}
}

\value{
\item{call}{The shared frailty model evaluated.}
\item{frequencies}{Numbers of patients, events and groups used to fit the model.}
\item{Nboot}{Number of bootstrap resamplings required.}
\item{Nbproblem}{Number of bootstrap resamplings in which the model did not converge. Those are excluded when computing the 
			standard-error and the percentile 95\% confidence interval.}
\item{ties}{A binary, indicating if the pairs tied on prediction were used to compute the concordance.}
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

\examples{

\dontrun{

## load data
data(readmission)

## a frailtypenal fit
fit.frailtypenal <- frailtyPenal(Surv(time,event)~dukes+cluster(id)+
                       charlson+chemo+strata(sex),data=readmission,
                       Frailty=TRUE,n.knots=10,kappa1=1,kappa2=1)

## a Cmeasures call
fit.Cmeasures <- Cmeasures(fit.frailtypenal)

## a short summary
print(fit.Cmeasures)

}

}

\keyword{methods}

