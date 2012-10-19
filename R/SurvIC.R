"SurvIC" <- function(lower,upper,event,tronc=FALSE,t0) {

if (missing(lower)) stop("Must have a lower time argument")
if (!is.numeric(lower)) stop("Lower time variable is not numeric")
if (any(is.na(lower))) stop("There is some NA in your lower time")

if (missing(upper)) stop("Must have an upper time argument")
if (!is.numeric(upper)) stop("Upper time variable is not numeric")
if (any(is.na(upper))) stop("There is some NA in your upper time")

if (missing(event)) stop("Must have an event argument")
if (!is.numeric(event)) stop("Event variable is not numeric")
if (any(is.na(event))) stop("There is some NA in your event argument")

if (length(lower)!=length(upper)) stop("Lower and upper time are different lengths")
if (any(lower>upper)) stop("Lower time has to be less than upper time")

if (length(lower==upper) != length(event==0)) warning("There may be an error in the right censored times")

if (!is.logical(tronc) || length(tronc)!=1) stop ("Tronc must be a simple logical")

if (tronc==FALSE) {
	if (!missing(t0)) stop("When tronc is FALSE, t0 argument must be deleted")
	ss <- cbind(lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "interval"
}
else {
	if (missing(t0)) stop("When tronc is TRUE, must have a troncature time t0")
	ss <- cbind(t0=t0,lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "intervaltronc"
}

if (is.R()) { class(ss) <- "SurvIC" } 
else { oldClass(ss) <- "SurvIC" }
return(ss)
}