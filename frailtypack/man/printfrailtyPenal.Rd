\name{print.frailtyPenal}
\alias{print.frailtyPenal}

\title{Print a Short Summary of parameter estimates of a shared gamma frailty model}
\description{
  Prints a short summary of 'frailtyPenal' object 
}
\usage{
print.frailtyPenal(x, digits = max(options()$digits - 4, 3), ...)
}
\arguments{
  \item{x}{the result of a call to the frailtyPenal function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\value{
 x, with the invisible flag set
}


\seealso{ \code{\link{summary.frailtyPenal}},  \code{\link{frailtyPenal}} }

\examples{
data(kidney)
mod<-frailtyPenal(Surv(time,status)~sex+age+cluster(id),
     n.knots=8,kappa1=10000,data=kidney)
print(mod)
}

\keyword{survival }
