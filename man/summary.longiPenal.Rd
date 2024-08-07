% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.longiPenal.R
\name{summary.longiPenal}
\alias{summary.longiPenal}
\alias{print.summary.longiPenal}
\title{Short summary of fixed covariates estimates of a joint model for
longitudinal data and a terminal event.}
\usage{
\method{summary}{longiPenal}(object, level = 0.95, len = 6, d = 2,
lab=c("coef","hr"), ...)
}
\arguments{
\item{object}{an object inheriting from \code{longiPenal} class}

\item{level}{significance level of confidence interval. Default is 95\%.}

\item{len}{the total field width for the terminal part. Default is 6.}

\item{d}{the desired number of digits after the decimal point. Default of 6
digits is used.}

\item{lab}{labels of printed results for the longitudinal outcome and the
terminal event respectively.}

\item{\dots}{other unused arguments.}
}
\value{
For the longitudinal outcome it prints the estimates of coefficients
of the fixed covariates with their standard error and p-values of the Wald
test. For the terminal event it prints HR and its confidence intervals for
each covariate. Confidence level is allowed (level argument).
}
\description{
This function returns coefficients estimates and their standard error with
p-values of the Wald test for the longitudinal outcome and hazard ratios
(HR) and their confidence intervals for the terminal event.  If a mediation analysis 
was performed (option \code{mediation} set to \code{TRUE} in \code{\link{longiPenal}}) 
this function displays estimations of the related quantities.
}
\examples{


\dontrun{
###--- Joint model for longitudinal data and a terminal event ---###

data(colorectal)
data(colorectalLongi)

# Survival data preparation - only terminal events 
colorectalSurv <- subset(colorectal, new.lesions == 0)

# Baseline hazard function approximated with splines
# Random effects as the link function

model.spli.RE <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS 
+ prev.resection, tumor.size ~  year * treatment + age + who.PS ,
colorectalSurv,	data.Longi = colorectalLongi, random = c("1", "year"),
id = "id", link = "Random-effects", left.censoring = -3.33, 
n.knots = 7, kappa = 2)

# Weibull baseline hazard function
# Current level of the biomarker as the link function

model.weib.CL <- longiPenal(Surv(time1, state) ~ age + treatment + who.PS
+ prev.resection, tumor.size ~  year * treatment + age + who.PS , 
colorectalSurv, data.Longi = colorectalLongi, random = c("1", "year"),
id = "id", link = "Current-level", left.censoring = -3.33, hazard = "Weibull")
	
summary(model.spli.RE)
summary(model.weib.CL)
}


}
\seealso{
\code{\link{longiPenal}}
}
\keyword{methods}
