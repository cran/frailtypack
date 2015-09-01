
\name{num.id}
\Rdversion{1.1}
\alias{num.id}
\title{
Identify individuals in Joint model for clustered data
}
\description{
This is a special function used in addition to the \code{cluster()} function in the context of survival joint models
for clustered data. This function identifies subject index. It is used on the right hand
side of a 'frailtyPenal' formula. Using \code{num.id()} in a formula implies that
a joint frailty model for clustered data is fitted (Rondeau et al. 2011).
}

\usage{
num.id(x)
}

\arguments{
\item{x}{A character or numeric variable which is supposed to indicate 
the variable identifying individuals}
}

\seealso{
 \code{\link{frailtyPenal}}
}

\references{
V. Rondeau, J.P. Pignon, S. Michiels (2011). A joint model for the dependance between clustered times to tumour progression and deaths: A meta-analysis of chemotherapy in head and neck cancer.
\emph{Statistical methods in medical research} \bold{897}, 1-19.

}

\examples{

\dontrun{

data(readmission)
#-- here is a generated cluster (31 clusters of 13 subjects)
readmission <- transform(readmission,group=id\%\%31+1)

###--- Joint Frailty model for clustered data ---###

#-- splines
joi.clus <- frailtyPenal(Surv(t.start,t.stop,event)~
cluster(group)+num.id(id)+dukes+charlson+sex+chemo+terminal(death),
formula.terminalEvent=~dukes+charlson+sex+chemo,
data=readmission,recurrentAG=TRUE,
n.knots=10,kappa=c(2.11e+08,9.53e+11))

}

}

\keyword{misc}
