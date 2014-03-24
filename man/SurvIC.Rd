\name{SurvIC}
\alias{SurvIC}

\title{
Create a survival object for interval censoring and possibly left truncated data
}
\description{
This is a function used in case of interval-censoring as a response variable in a model formula only for Cox proportional hazard or shared frailty model. Sometimes, an unobserved event might occur in a time interval [L,U]. RecurrentAG argument gets invalid with the use of SurvIC. Time dependant variables are not allowed.
}
\usage{
SurvIC(t0, lower, upper, event)
}

\arguments{
  \item{lower}{
Starting time of the interval for interval-censored data. Time of right-censoring instead.
}
  \item{upper}{
Ending time of the interval for interval-censored data. For right-censored data, lower and upper time must be equal (for numerical reason).
}
  \item{event}{
Status indicator 0=right-censored, 1=interval-censored
}
  \item{t0}{
Truncation time for left truncated data only. To be ignored otherwise.
}
}

\details{
Typical usages are
\code{SurvIC(lower,upper,event)} or \code{SurvIC(t0,lower,upper,event)}
}

\examples{

\dontrun{

data(bcos)
bcos$event <- ifelse(bcos$left!=bcos$right,1,0)

###---  Cox proportional hazard model with interval censoring ---###

cox.ic <- frailtyPenal(SurvIC(left,right,event)~treatment,
data=bcos,n.knots=8,kappa1=10000)

###---  Shared model with interval censoring ---###

bcos$group <- c(rep(1:20,4),1:14)

sha.ic <- frailtyPenal(SurvIC(left,right,event)~cluster(group)+
treatment,data=bcos,n.knots=8,kappa1=10000)

}

}
