
\name{subcluster}
\Rdversion{1.1}
\alias{subcluster}
\title{
Identify subclusters
}
\description{
This is a special function used in the context of survival nested models.  It
identifies correlated groups of observations within other groups defined by using 
'cluster' function from 'survival' package, and is used on the right hand
side of 'frailtyPenal' formula for fitting a nested model.  Using \code{subcluster()} in a formula implies that
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
 \code{\link{frailtyPenal}}
}

\examples{

\dontrun{

data(dataNested)
modClu <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+
subcluster(subgroup)+cov1+cov2,data=dataNested,
n.knots=8,kappa1=50000,kappa2=50000,hazard="Splines")

print(modClu)

}

}

\keyword{misc}
