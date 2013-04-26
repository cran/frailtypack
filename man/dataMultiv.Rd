\name{dataMultiv}
\alias{dataMultiv}
\docType{data}
\title{
Simulated data with two types of recurrent events and a dependant terminal event}

\description{
This contains simulated data with two different types of recurrent event times and a terminal event time from a multivariate frailty models with two gaussian random effects. This data set is useful to illustrate how to fit a multivariate frailty model.
This database combines different observations for each subject (sample size=250, number of lines=1237 so in average four recurrent events per subject). The right-censored case is treated only and the calendar timescale representation is used. Exponential death times and recurring times are generated using a gaussian  multivariate frailty model. The frailty terms are iid \eqn{N(0,0.5)} and correlated with \eqn{\rho=0.5} and \eqn{\alpha1=0.6} and \eqn{\alpha2=-0.3}, the coefficient associated with the frailty term in the death time hazard function, is set at 0. The right-censoring variable is set at 5 and the covariate is generated from a Bernoulli distribution (p=0.5).  
}

\usage{data(dataMultiv)}
\format{
  A data frame with 1237 observations on the following 9 variables.
  \describe{
    \item{\code{t0}}{start of interval (0 or previous recurrence time).}
    \item{\code{t1}}{end of interval (recurrent, censoring or terminal event time).}
    \item{\code{id}}{identification variable (ex: subject number).}
    \item{\code{deltaEvent1}}{censoring status for the first type of recurrent events (1:event, 0: censored or dead.}
    \item{\code{deltaEvent2}}{censoring status for the second type of recurrent events (1:event, 0: censored or dead.}
    \item{\code{deltadc}}{1:dead status or 0:aline.}
    \item{\code{var1}}{dichotomous covariate (0,1).}
    \item{\code{var2}}{dichotomous covariate (0,1).}
    \item{\code{var3}}{dichotomous covariate (0,1).}
  }
}

\keyword{datasets}
