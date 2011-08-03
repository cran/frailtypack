
\name{dataJoint}
\docType{data}
\alias{dataJoint}

\title{Simulated data with recurrent events and informative terminal event}
\usage{data(dataJoint)}

\description{
This contains simulated data with recurrent events and time until last follow-up that is related to the recurrent process. This data set is useful to illustrate how to fit a joint model where the censoring mechanism is informative about the recurring process.
This database combines different observations for each subject (sample size=500, number of lines=1053 so in average more than two recurrent events per subject). The right-censored case is treated only and the calendar timescale representation is used with delayed entry. Exponential death times and recurring times are generated using a gamma frailty model. The common frailty term is iid \eqn{\Gamma(2,2)} and \eqn{\alpha}, the coefficient associated with the frailty term in the death time hazard function, is set at 0. The right-censoring variable is set at 0.8 and the two covariates are generated from a Bernoulli distribution (p=0.5).  

}

\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{identification variable (ex: subject number) }
    \item{time.entry}{start of interval (0 or previous recurrence time)}
    \item{time.end}{end of interval (recurrence or censoring time)}
    \item{status}{censoring status for recurrent events (1:event, 0: censored or dead}
    \item{status.terminal}{censoring status for terminal events (0:alive, 1:death)}
    \item{var1}{dichotomous covariate (0,1)}
    \item{var2}{dichotomous covariate (0,1)}
  }
}
\seealso{

\code{\link{frailtyPenal for Joint frailty models}}

}
 

\source{

V. Rondeau, S. Mathoulin-Pellissier, H. Jacqmin-Gadda, V. Brouste, P. Soubeyran (2007). Joint frailty models for recurring events and death using maximum penalized likelihood estimation:application on cancer events.
\emph{Biostatistics}, \bold{8},4, 708-721.

}

\keyword{datasets}
