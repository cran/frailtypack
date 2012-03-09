
\name{event2}
\Rdversion{1.1}
\alias{event2}
\title{
Identify event2 indicator 
}
\description{
This is a special function used in the context of recurrent event models with event2 event
(e.g., censoring variable related to recurrent events). It contains the status indicator, normally 0=alive, 1=dead, and is used on the right hand side of a formula  of a 'frailtyPenal' object.  Using \code{event2()} in a formula implies that a joint frailty model for recurrent events and event2 events is fitted.
}

\usage{
event2(x)

}
\arguments{
\item{x}{
A numeric variable but should be a Boolean which equals 1 if the subject is dead and 0 if he is alive or censored, as a death indicator.
}}
\value{
\item{x}{
a death indicator}

}


\seealso{ \code{\link{frailtyPenal for Joint frailty models}}, \code{\link{plot.jointPenal}}, \code{\link{print.jointPenal}}, \code{\link{summary.jointPenal}}
}


\examples{


\dontrun{

  data(dataJoint)

  modJoint<-frailtyPenal(Surv(time.entry,time.end,status)~cluster(id)+var1+var2
                         +terminal(status.terminal),
                         formula.terminalEvent=~var1,
                         data=dataJoint,n.knots=7,Frailty=TRUE,
                         kappa1=1, kappa2=1, joint=TRUE, recurrentAG=TRUE)
    
  print(modJoint)
}


}

\keyword{misc}






