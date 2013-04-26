
\name{event2}
\Rdversion{1.1}
\alias{event2}
\title{
Identify event2 indicator 
}
\description{
This is a special function used in the context of multivariate frailty model with two types of recurrent events and a terminal event
(e.g., censoring variable related to both recurrent events). It contains the indicator of the recurrent event of type 2, normally 0=no event, 1=event, and is used on the right hand side of a formula  of a 'multivePenal' object.  Using \code{event2()} in a formula implies that a multivariate frailty model for two types of recurrent events and a terminal event is fitted.
}

\usage{
event2(x)

}
\arguments{
\item{x}{
A numeric variable but should be a Boolean which equals 1 if the subject has experienced an event of type 2 and 0 if not.
}}
\value{
\item{x}{
a event of type 2 indicator}

}


\seealso{ \code{\link{multivePenal}} }


\examples{

\dontrun{

data(dataMultiv)

ModMultive <- multivePenal(Surv(t0,t1,deltaEvent1)~
cluster(id)+var1+var2+terminal(deltadc)+event2(deltaEvent2),
formula.terminalEvent=~var1,formula2=~var1+var2+var3,
data=dataMultiv,Frailty=TRUE,recurrentAG=TRUE,
nb.int=c(3,3,3),hazard="Piecewise-per")

print(ModMultive)

}

}

\keyword{misc}






