
\name{subcluster}
\Rdversion{1.1}
\alias{subcluster}
\title{
Identify subclusters
}
\description{
This is a special function used in the context of survival nested models.  It
identifies correlated groups of observations within other groups defined
 by using 
'cluster' function from 'survival' package, and is used on the right hand
side of a formula of a 'frailtyPenal' object for fitting a nested model.  Using \code{subcluster()} in a formula implies that
a nested frailty model is estimated. 
}

\usage{
subcluster(x) 
}

\arguments{
\item{x}{A character, factor, or numeric variable which is supposed to indicate 
the variable subgroup}
}

\value{
\item{x}{A variable identified as a subcluster }
}

\seealso{
\code{\link{frailtyPenal for Nested frailty models}},  \code{\link{print.nestedPenal}},
\code{\link{plot.nestedPenal}},  \code{\link{summary.nestedPenal}}
}

\examples{

\dontrun{

  data(dataNested)
  modClu<-frailtyPenal(Surv(t1,t2,event)~cluster(group)+subcluster(subgroup)+
            cov1+cov2,Frailty=TRUE,data=dataNested,n.knots=8,kappa1=50000,kappa2=50000)

  print(modClu)

# It takes aound 24 minutes to converge (depends on the processor) #
}

}

\keyword{misc}
