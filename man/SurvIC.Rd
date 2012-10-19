\name{SurvIC}
\alias{SurvIC}

\title{
Create a survival object for interval censoring and possibly left truncated data
}
\description{
This is a function used in case of interval-censoring as a response variable in a model formula only for Cox proportional hazard or shared frailty model. Sometimes, an unobserved event might occur in a time interval [L,U]. RecurrentAG argument gets invalid with the use of SurvIC. Time dependant variables are not allowed.
}
\usage{
SurvIC(lower, upper, event, tronc = FALSE, t0)
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
  \item{tronc}{
Truncation indicator, logical value. There is a possiblity to add a truncature time, if TRUE then t0 is needed. The default is FALSE.
}
  \item{t0}{
Truncation time. Must be deleted if no left truncation.
}
}

\examples{

data(diabetes)
# here is created the needed variable event
diabetes <- transform(diabetes,event=ifelse(left!=right,1,0))
# here is created a hypothetical groupement
diabetes$groupe <- rep(1:146,c(rep(5,145),6))

###  Cox proportional hazard model with interval censoring  ###

frailtyPenal(SurvIC(left,right,event)~gender,
             data=diabetes,Frailty=FALSE,n.knots=8,cross.validation=FALSE,
             kappa1=3309,intcens=TRUE)

###  Shared model with interval censoring  ###

frailtyPenal(SurvIC(left,right,event)~cluster(groupe)+gender,
             data=diabetes,Frailty=TRUE,n.knots=8,cross.validation=FALSE,
             kappa1=3309,intcens=TRUE)

}