
\name{plot.prediction}
\alias{plot.predFrailty}
\alias{plot.predJoint}

\title{Plot Method for predictions.}

\description{
  Probabilites of events and deaths. Confidence bands are allowed. 
}

\usage{
  \method{plot}{predFrailty}(x, conf.bands=FALSE, pos.legend="topleft", cex.legend=0.7, ...)
}

\arguments{
  \item{x}{An object from the 'prediction' function}
  \item{conf.bands}{Logical value. Determines whether confidence bands will be plotted. The default is FALSE.}
  \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the list
                    '"bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and '"center"'. The default is '"topright"'}
  \item{cex.legend}{size of the legend. Default is 0.7}
  \item{\dots}{Other unused arguments}
}

\value{
  Print as many plots as there are individuals with their probabilities to die.
}

\examples{

\dontrun{

data(readmission)

##-- predictions on a Cox proportionnal hazard model --##
cox <- frailtyPenal(Surv(time,event)~sex+dukes,
n.knots=10,kappa1=10000,data=readmission,Frailty=FALSE)

#-- construction of the dataframe for predictions
datapred <- data.frame(sex=0,dukes=0)
datapred$sex <- as.factor(datapred$sex)
levels(datapred$sex)<- c(1,2)
datapred$dukes <- as.factor(datapred$dukes)
levels(datapred$dukes)<- c(1,2,3)
datapred[1,] <- c(1,2)
datapred[2,] <- c(2,3)

#-- prediction of the event for two patients, between time 100 and time 2000
pred.cox <- prediction(cox,datapred,100,50,2000)
plot(pred.cox)

}

}

\keyword{file}
