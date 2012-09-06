
\name{diabetes}
\docType{data}
\alias{diabetes}

\title{Interval-censored data for time from onset of diabetes to the onset of diabetic nephronpathy}
\usage{data(diabetes)}

\description{
 The diabetes dataframe is based on a study conducted at the Steno Memorial Hospital in Denmark from 1933-1984. 
The data involve 731 patients (454 males and 277 females) and have 595 exact, 136 interval-censored and no right-censored observations.
}

\format{
  A data frame with 731 observations on the following 3 variables:
  \describe{
    \item{left}{left end point of a censoring interval}
    \item{right}{right end point of a censoring interval}
    \item{gender}{gender indicator (1 for male and 0 for female)}
  }
}

\source{
 Kindly provided by Professor P.K. Andersen.
}

\keyword{datasets}
