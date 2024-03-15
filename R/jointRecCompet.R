#' Competing Joint Frailty Model: A single type of recurrent event and two
#' terminal events.
#'
#' @description Fit a joint competing frailty model for a single recurrent event
#' and two terminal events defined as,
#' 
#' \deqn{\text {Recurrent event:} \quad r_{i j}\left(t \mid w_i, 
#' \boldsymbol{X}_{r, i j}\right)=r_0(t) \exp \left(\boldsymbol{X}_{r, i j} \boldsymbol{\beta}_r+w_i\right)}
#' \deqn{\text{First terminal event:} \quad \lambda_{1, i}\left(t \mid w_i, 
#' \boldsymbol{X}_{1 i}\right)=\lambda_{1, 0}(t) \exp \left(\boldsymbol{X}_{1, i} \boldsymbol{\beta}_1+\alpha_1 w_i\right)}
#' \deqn{\text{Second terminal event:} \quad \lambda_{2, i}\left(t \mid w_i, 
#' \boldsymbol{X}_{2, i}\right)=\lambda_{2, 0}(t) \exp \left(\boldsymbol{X}_{2, i} \boldsymbol{\beta}_2+\alpha_2 w_i\right).}
#' 
#' where \eqn{\omega_i \sim \mathcal{N}(0,\theta)} is the frailty term and \eqn{\boldsymbol{X}_{r, i j},\boldsymbol{X}_{1, i}} and
#' \eqn{\boldsymbol{X}_{2, i}} are vectors of baseline covariates (possibly the same). The parameters \eqn{\alpha_1} and \eqn{\alpha_2}
#' are power parameters.
#' 
#' @aliases jointRecCompet 
#'
#' @usage
#'  jointRecCompet(formula,
#'          formula.terminalEvent = NULL,
#'          formula.terminalEvent2 = NULL,
#'          data,
#'          initialize = TRUE,
#'          recurrentAG = FALSE,
#'          maxit = 350,
#'          hazard = "Weibull",
#'          n.knots=7,
#'          kappa = rep(10, 3),
#'          crossVal=FALSE,
#'          constraint.frailty = "squared",
#'          GHpoints = 32,
#'          tolerance = rep(10^-3, 3),
#'          init.hazard = NULL,
#'          init.Sigma = 0.5,
#'          init.Alpha1 = 0.1,
#'          init.Alpha2 = -0.1,
#'          init.B = NULL)
#'
#' @param formula a formula object, with the response for the first recurrent
#' event on the left of a \eqn{\sim} operator, and the terms on the right.
#' The response must be in the format Surv(t0, t1, recurrentevent) for a calendar-time
#' specification where \code{t0} is the start time for an at-risk period for the recurrent event,
#' \code{t1} is the end time for an at-risk period for the recurrent event, and
#' \code{recurrentevent} is a numeric indicator for whether an event was observed (1)
#' or was censored (2).
#' In a gap-time setting, an object of the format Surv(t, recurrentevent) should be used instead.
#' Note that to not be confused with a left-truncation setting, when using a calendar-time specification 
#' argument \code{recurrentAG} should be set to \code{TRUE}. 
#' 
#' @param formula.terminalEvent, a formula object,
#' empty on the left of a \eqn{\sim} operator,
#' and the terms on the right. Leave the formula at
#' the default value (NULL) for a model with no variables.
#'
#' @param formula.terminalEvent2, a formula object,
#' empty on the left of a \eqn{\sim} operator,
#' and the terms on the right.Leave the formula at
#' the default value (NULL) for a model with no variables.
#'
#' @param data a 'data.frame' with the variables used in 'formula',
#' 'formula.terminalEvent', and 'formula.terminalEvent2'.
#'
#' @param initialize Logical value to internally initialize regression coefficients and
#' baseline hazard functions parameters using simpler models from frailtypack.
#' When initialization is requested, the
#' program first fits two joint frailty models for the recurrent
#' events and each terminal event.
#' When FALSE, parameters are initialized via the arguments
#' init.hazard, init.Sigma, init.Alpha1, init.Alpha2, init.B.
#'
#' @param recurrentAG Logical value. Is Andersen-Gill model fitted? 
#' If so indicates that recurrent event times with the counting 
#' process approach of Andersen and Gill is used. This formulation can be 
#' used for dealing with time-dependent covariates. The default is FALSE.
#'
#' @param maxit maximum number of iterations for the Marquardt algorithm.
#' Default is 350.
#'
#' @param hazard Type of hazard functions. Available options are \code{"Weibull"} 
#' for parametric Weibull function, \code{"Splines"} for semiparametric
#' hazard functions using equidistant intervals or \code{"Splines-per"} for
#' percentile intervals. Default is \code{"Weibull"}.
#'
#' @param n.knots In the case of splines hazard functions, number of knots to be used 
#' in the splines basis. This number should be between 4 and 20. Default is 7. 
#' 
#' @param kappa In the case of splines hazard functions, a vector of size 3 containing
#' the values of the smoothing parameters to be used for each baseline hazard function. 
#' Default value is 10 for each function. 
#' 
#' @param crossVal In the case of splines hazard functions, indicates how 
#' the smoothing parameters are chosen. If set to "TRUE" then those parameters
#' are chosen automatically using cross-validation on reduced models for each 
#' baseline hazard function. If set to "FALSE" then the parameters are those provided 
#' by the argument \code{kappa}. If set to "TRUE" then the argument \code{kappa} is 
#' ignored. Default is "TRUE". 
#' 
#' 
#' @param constraint.frailty Type of positivity constraint used for the variance of
#' of the random effect in the likelihood. Possible values are 'squared' or 'exponential'. 
#' Default is 'squared'. See Details. 
#' 
#' @param GHpoints Integer. Number of nodes for Gauss-Hermite integration
#' to marginalize random effects/frailties. Default is 32.
#'
#' @param tolerance Numeric, length 3. Optimizer's tolerance for (1) successive change
#' in parameter values, (2) log likelihood, and (3) score, respectively.
#'
#' @param init.hazard Numeric. Initialization values for hazard parameters.
#' If a weibull model is used, the order is:
#' shapeR, scaleR, shapeTerminal1, scaleTerminal1, shapeTerminal2, scaleTerminal2.
#'
#' @param init.Sigma Numeric,. Initialization value for the standard deviation of the
#' normally-distributed random effects.
#'
#' @param init.Alpha1 Numeric. Initialization value for the parameter alpha that
#' links the hazard function of the recurrent event to the first terminal event.
#'
#' @param init.Alpha2 Numeric. Initialization value for the parameter alpha that
#' links the hazard function of the recurrent event to the second terminal event.
#'
#' @param init.B Numeric vector of the same length and order 
#' as the three covariate vectors for the recurrent, terminal1,
#' and terminal2 events (in that order).
#'
#' @details
#' Right-censored data are allowed.
#' Left-truncated data and stratified analysis are not possible. 
#' Prediction options are not yet available. 
#' The \code{constraint.frailty} argument defines the positivity constraint 
#' used for the frailty variance in the likelihood. By default it uses the square 
#' so that the absolute value of the parameter is the standard deviation of the frailty
#' (i.e \eqn{\theta^2 = \beta^2}). 
#' The other parametrization uses the square of the exponential for the variance
#' so that the parameter is the logarithm of the standard deviation (\eqn{\theta^2 = (\exp{\beta})^2}).
#' For others parameters in the model needing a positivity constraint (parameters related to the
#' baseline hazard functions), the parametrization used is the exponential squared.
#' 
#'
#' @return Parameters estimates of a competing joint frailty model, more
#' generally a 'jointRecCompet' object. Methods defined for 'jointRecCompet' objects
#' are provided for print, plot and summary. The following components are
#' included in a 'jointRecCompet' object.
#' 
#' \item{summary.table}{A table describing the estimate, standard error,
#' confidence interval, and pvalues for each of the parameters in the model.}
#' 
#' \item{controls}{A vector of named control parameters}
#' 
#' \item{k0}{For splines baseline hazard functions, vector of penalization terms.}
#' 
#' \item{noVarEvent}{A vector containing for each event type if there is no covariate used
#' in the model.}
#' 
#' \item{np}{Total number of parameters}
#' 
#' \item{b}{Vector containing the estimated coefficients of the model before any positivity constraint.
#' The values are in order: the coefficients associated with the baseline hazard functions 
#' (either the splines or the shape and scale parameters for Weibull hazard), 
#' the random effect variance, the coefficients of the frailty (\eqn{\alpha_1} and \eqn{\alpha_2}) 
#' and the regression coefficients.}
#'
#' \item{H_hessOut}{Covariance matrix of the estimated parameters}
#' 
#' \item{HIHOut}{Covariance matrix of the estimated parameters for the penalized likelihood
#' in the case of Splines baseline hazard functions.}
#' 
#' \item{LCV}{The approximated likelihood cross-validation criterion in the spline case}
#'
#' \item{critCV}{Convergence criteria}
#' 
#' \item{x1}{Vector of times for which the hazard function of
#' the recurrent event is estimated.
#' By default seq(0,max(time),length=99),
#' where time is the vector of survival times.}
#'
#' \item{lam1}{Matrix of hazard estimates and confidence bands for the recurrent
#' event.}
#'
#' \item{xSu1}{Vector of times for the survival function of
#' the recurrent event.}
#'
#' \item{surv1}{Matrix of baseline survival
#' estimates and confidence bands for recurrent event.}
#' 
#' \item{x2}{Vector of times for the first terminal event (see x1 value).}
#'
#' \item{lam2}{Matrix of hazard estimates and confidence bands for the
#' first terminal event.}
#' 
#' \item{xSu2}{Vector of times for the survival function of the first terminal event.}
#'
#' \item{surv2}{Vector of the survival function of the first terminal event evaluated at xSu2.}
#'
#' \item{x3}{Vector of times for the second terminal event (see x1 value).}
#'
#' \item{lam3}{Matrix of hazard estimates and confidence bands for the
#' second terminal event.}
#'
#' \item{xSu3}{Vector of times for the survival function of the second terminal event.}
#'
#' \item{surv3}{Vector of the survival function of the second terminal event evaluated at xSu3.}
#' 
#' \item{ni}{Number of iterations needed to converge.}
#'
#' \item{constraintfrailty}{Positivity constraint used for the variance of the random effect} 
#' 
#' \item{ziOut1}{In the spline case, vector of knots used in the spline basis for the recurrent event}
#' 
#' \item{ziOutdc}{In the spline case, vector of knots used in the spline basis for the terminal events}
#' 
#' \item{ghnodes}{Nodes used for the Gauss-Hermite quadrature.}
#' 
#' \item{ghweights}{Weights used for the Gauss-Hermite quadrature.}
#' 
#' \item{tolerance}{Numeric, length 3. Optimizer's tolerance for (1) successive change
#' in parameter values, (2) log likelihood, and (3) score, respectively.}
#' 
#' \item{call}{Call of the function.}
#' 
#' \item{loglikPenal}{Estimated penalized log-likelihood in the spline case}
#' 
#' \item{logLik}{Estimated log-likelihood in the Weibull case}
#' 
#' \item{AIC}{For the Weibull case, Akaike Information criterion}
#' 
#' \item{n}{Total number of subjects}
#' 
#' \item{nevts}{Number of events for each event type.}
#'  
#' @seealso
#' \code{\link{terminal}}
#'
#' @export
#' @importFrom stats drop.terms
#' @examples 
#' \donttest{
#' set.seed(1)
#' data=simulatejointRecCompet(n=500,
#'			par0=c(shapeR = 1.5, scaleR = 10, 
#'			shapeM = 1.75, scaleM = 16, shapeD = 1.75, scaleD = 16, sigma = 0.5, 
#'			alphaM = 1, alphaD = 1, betaR = -0.5, betaM = -0.5, betaD = 0) )
#' mod <-jointRecCompet(formula = Surv(tstart, tstop, event)~cluster(id)+treatment+
#'                     terminal(terminal1)+terminal2(terminal2),
#'                     formula.terminalEvent = ~treatment,
#'                     formula.terminalEvent2 = ~treatment,
#'                     data = data,
#'                     recurrentAG = TRUE,
#'                     initialize = TRUE,
#' 					   n.knots=7,
#' 					   crossVal=TRUE,
#'                     hazard = "Splines",
#'                     maxit = 350)
#' 
#' #This example uses an extract of 500 patients of the REDUCE trial
#' data(reduce)
#' mod_reduce <-jointRecCompet(formula = Surv(t.start,t.stop, del)~cluster(id)+
#' 					   treatment+terminal(death)+terminal2(discharge),
#'                     formula.terminalEvent = ~treatment,
#'                     formula.terminalEvent2 = ~treatment,
#'                     data = reduce,
#'                     initialize = TRUE,
#'                     recurrentAG = TRUE,
#'                     hazard = "Weibull",
#'                     constraint.frailty = "exponential",
#'                     maxit = 350)
#' print(mod_reduce)
#' }

"jointRecCompet" <-
  function(formula,
           formula.terminalEvent = NULL,
           formula.terminalEvent2 = NULL,
           data,
           initialize = TRUE,
           recurrentAG = FALSE,
  		   maxit = 350,
           hazard = "Weibull",
		   n.knots=7,
		   kappa = rep(10, 3),
		   crossVal=FALSE,
		   constraint.frailty="squared",
           GHpoints = 32,
           tolerance = rep(10^-3, 3),
           init.hazard = NULL,
           init.Sigma = 0.5,
           init.Alpha1 = 0.1,
           init.Alpha2 = -0.1,
           init.B = NULL)
{
#############################################################
#############################################################
# Outline of competingPenal.R
# (1) Verify/extract event indicators, times, groups
# (2) Configure Hazards
# (3) Configure Model Matrices
# (4) Configure Parameters
# (5) Define Gauss-Hermite Nodes and Weights
# (6) Fill Starting Parameter Vector with User-Defined Values
#		OR Initialize Models
# (7) Check Dimensions of all variables being sent to Fortran
#       for debugging
# (8) Send to Fortran for optimization
# (9) Format Model Summary Table
# (10) Format Initialization Model Summary Tables
#		(if initialization == TRUE)
#############################################################
#############################################################

#n.knots <- 4 # remove when splines models are available
#kappa <- rep(5, 3)
#crossVal <- 0

jointGeneral <- FALSE
if(!inherits(formula,"formula")) stop("The argument formula must be a formula")

if (nrow(data) == 0)
  stop("No (non-missing) observations")
#############################################################
if(jointGeneral) stop("Model with correlated random effects not yet implemented.")
#############################################################
# (1) Verify/extract event indicators, times, groups

### (1a) Verify Recurrent Event 1 / Times
# NN =  names of time gap and recurrent event 1
Y1 <- get_all_vars(update(formula, "~1"), data)

if (ncol(Y1) == 3){
	if(!recurrentAG) stop("When using a calendar-time format, 'recurrentAG' must be set to 'TRUE'. ")
	cat('Using a calendar-time format','\n')
	gapTimes=FALSE
}else if(ncol(Y1)==2){
	cat('Using a gap-time format','\n')
	gapTimes=TRUE
}else{
	stop('Error in the specification of the survival object.')
}
if(!(constraint.frailty %in% c('squared','exponential'))){
	stop("'constraint.frailty' must be 'squared' or 'exponential'.")
}else{
	constraint.frailty<-switch(constraint.frailty,
  	     "squared" = 0,
		 "exponential"= 1)
}
NN <- colnames(Y1)

if(!gapTimes){
	TSTART <- NN[1]
	TSTOP <- NN[2]
	EVENT1 <- NN[3]
	tt10 <- Y1[, 1]
	tt11 <- Y1[, 2]
	event1 <- Y1[, 3]
}else{
	TSTART <- 'tstart' 
	TSTOP <- NN[1]
	EVENT1 <- NN[2]
	tt11 <- Y1[, 1]
	event1 <- Y1[, 2]
	tt10 <- rep(0,length(tt11))

}


call=match.call()

### (1b) Verify recurrent event 2 (Not yet implemented)
event2.ind <- 0
data.Event2 <- NULL
formula.Event2 <- NULL
event2 <- 0
tt0meta0 <- 0
tt1meta0 <- 0
kappa0 <- 0
 


### (1c) Verify Terminal Event 1
TT <-
  survival::untangle.specials(terms(formula, c("terminal")), "terminal", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
TERMINAL1 <- substr(TT, start = start + 1, stop = stop - 1)
if (length(TERMINAL1) == 0) {
  stop("A term for a terminal event must be included in the formula.")
}
if (!all(data[[TERMINAL1]] %in% c(1, 0))) {
  stop("terminal must contain a variable coded 0-1 and a non-factor variable")
}

### (1d) Verify Terminal Event 2
TT <-survival::untangle.specials(terms(formula, c("terminal2")), "terminal2", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
TERMINAL2 <- substr(TT, start = start + 1, stop = stop - 1)
if (length(TERMINAL2) == 0) {
  terminal2.ind <- 0
} else{
  if (!all(data[[TERMINAL2]] %in% c(1, 0))) {
  	stop(
  		"terminal must contain a variable coded 0-1 and a non-factor variable"
  	)
  }
  terminal2.ind <- 1
}

### (1e) Verify Cluster Variable
# name of cluster variable
TT <-
  survival::untangle.specials(terms(formula, c("cluster")), "cluster", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
CLUSTER <- substr(TT, start = start + 1, stop = stop - 1)

if (length(data[[CLUSTER]])) {
  uni.cluster <- unique(data[[CLUSTER]])
} else{
  stop("grouping variable is needed")
}
if (length(uni.cluster) == 1) {
  stop("grouping variable must have more than 1 level")
}
if (event2.ind == 1){
	if(CLUSTER %in% colnames(data.Event2)) {
		stop("grouping variable must be present in data.Event2")
	}
	if (!all(uni.cluster %in% data.Event2[[CLUSTER]])) {
		stop("all groups must be represented in data.Event2")
	}
}
##########################################################################
##########################################################################
# (2) Configure Hazards
if (!(hazard %in% c("Weibull","Splines","Splines-per"))){
  stop("Argument 'hazard' must be one of the following: 'Weibull', 'Splines' or 'Splines-per'.")
}

haztemp <- hazard

# Specify whether spline points are regularly spaced or at percentiles
if(hazard %in% c("Splines","Splines-per")){
	equidistant <- ifelse(hazard=="Splines",1,0)
}else{
  	### Weibull
  	equidistant <- 1
}
# typeof is the numerical indicator for hazard
typeof <- switch(hazard,
  	     "Splines" = 0,
		 "Splines-per"=0,
  	     "Weibull" = 2)


#### Configure Splines Hazard (knots, penalty), Not yet implemented
 if (typeof == 0) {
   #crossVal <- 1 # always do cross validation for splines
   #if (missing(kappa) & crossVal==F)
   #	stop("smoothing parameter (kappa1) is required")
   #if (missing(n.knots))
   #	stop("number of knots are required")

   if(!inherits(kappa,"numeric")) stop("The argument kappa must be a numeric")
   	
   if(!inherits(n.knots,"integer")) n.knots <- as.integer(n.knots)

   if (length(n.knots) != 1) {
   	stop("length of knots must be 1.")
   }
   if (length(kappa) != 2 + event2.ind + terminal2.ind) {
   	stop(
   		"length of kappa must be equal to the number of formulas
   		for different event types (3 or 4) in the order
   		recurrent1, recurrent2, terminal1, terminal2."
   	)
   }
   if (event2.ind == 1 & terminal2.ind == 0) {
   	kappa = c(kappa0, 0)
   }
   if (event2.ind == 0 & terminal2.ind == 1) {
   	kappa = c(kappa[1:2], 0, kappa[3])
   }
   n.knots[n.knots < 4 & n.knots != 0] <- 4
   n.knots[n.knots > 20] <- 20
 }else if(typeof == 2){
   # if the hazard is weibull
   if (!(missing(n.knots)) || !(missing(kappa))) {
   	warning("When parametric hazard is not 'Splines'
   		'kappa' and 'n.knots' arguments are ignored.")
   }
   n.knots <- 0
   kappa <- rep(0, 4)
   crossVal <- 0
 }

# End Hazard Configuration
#########################################################################
#########################################################################
# (3) Configure Model Matrices

# noVarEvent indicates whether there are no explanatory variables for the
# recurrent1, terminal1, recurrent2, and terminal2 events (in that order)
noVarEvent = c(0,0,1,0)

# delete specials from the formula to get just covariates
# on right side.
specials = c("strata", "cluster", "terminal", "event2", "terminal2")
Terms = terms(formula, specials = specials)

# Check if all terms on right side of Recurrent event formula are specials
if(length(unlist(attr(Terms, "specials"))) == (length(unlist(attr(Terms, "variables")))-2)){
	noVarEvent[1] <- 1
}

# Check for terminal1 formula
if(is.null(formula.terminalEvent)){
	noVarEvent[2] <- 1
}

# Check for terminal2 formula
if(is.null(formula.terminalEvent2)){
	noVarEvent[4] <- 1
}

# Recurrent Event 1 Model Matrix
if(noVarEvent[1] == 0){
modelmatrix1 =
  model.matrix(update(
  	drop.terms(
  		termobj = terms(formula),
  		unlist(attr(Terms, "specials")) - 1,
  		keep.response = TRUE
  	),
  	~ . - 1
  ),
  data)
}else{
	  modelmatrix1 = matrix(0)
}

# need to get densely-ranked ids
group1 <- as.numeric(factor(data[[CLUSTER]]))

# Compute Event Counts
# nevents1 <- tapply(event1, group1, sum)

# Recurrent Event 2 Model Matrix
if (event2.ind == 1) {
  group2 <- as.numeric(factor(data.Event2[[CLUSTER]]))
  # Compute Event Counts
  #nevents2 <- tapply(event2, group2, sum)
  Terms2 = terms(formula.Event2, specials = specials)
  if (!is.null(unlist(attr(Terms2, "specials")))) {
  	modelmatrix3 =
  		model.matrix(update(
  			drop.terms(
  				termobj = terms(formula.Event2),
  				unlist(attr(
  					Terms2, "specials"
  				)) - 1,
  				keep.response = TRUE
  			),
  			~ . - 1
  		),
  		data.Event2)
  } else{
  	modelmatrix3 =
  		model.matrix(update(formula.Event2, ~ . - 1),
  			 data.Event2)
  }
} else{
  modelmatrix3 = matrix(0)
  group2 = 0
}

# Terminal Event 1 Model Matrix
data.terminal <- do.call(what = "rbind",
  	 lapply(split(x = data, f = data[[CLUSTER]]),
  	        function(df) {
  	        	subset(df, df[[TSTOP]] == max(df[[TSTOP]]))
  	        }))
groupdc <- as.numeric(factor(data.terminal[[CLUSTER]]))
tt1dc <- data.terminal[[TSTOP]]
terminal1 <- data.terminal[[TERMINAL1]]

if(noVarEvent[2]==0){
modelmatrix2 = model.matrix(update(formula.terminalEvent, ~ . - 1), data.terminal)
}

# Terminal 2 Model Matrix
if((terminal2.ind == 1)) {
	if(noVarEvent[4] == 0){
		  modelmatrix4 = model.matrix(update(formula.terminalEvent2, ~ . - 1),
  		    data.terminal)
	}else{
		  modelmatrix4 = matrix(0)
	}
    terminal2 <- data.terminal[[TERMINAL2]]
} else{
    terminal2 <- data.terminal[[TERMINAL1]]*0
    modelmatrix4 = matrix(0)
}
#########################################################################
#########################################################################
# (4) Configure Parameters

### Total number of parameters
nvar = ncol(modelmatrix1) * (1-noVarEvent[1]) +
  ncol(modelmatrix2) * (1-noVarEvent[2])  +
  ncol(modelmatrix3) * event2.ind * (1-noVarEvent[3])+
  ncol(modelmatrix4) * terminal2.ind * (1-noVarEvent[4])

nbvar = c(
	ncol(modelmatrix1) * (1-noVarEvent[1]) ,
	ncol(modelmatrix2) * (1-noVarEvent[2]) ,
	ncol(modelmatrix3)* event2.ind * (1-noVarEvent[3]),
	ncol(modelmatrix4)* terminal2.ind * (1-noVarEvent[4])
)
# Total number of parameters
# This will need adjustment before incorporating a second recurrent event
if(typeof==0 ){ # splines, single random effect
	np = ((2 + n.knots) * (2 + event2.ind + terminal2.ind) + nvar + 3 + 2*jointGeneral)
}else if(typeof == 2){ # weibull, single random effect
	np = (2 * (2 + event2.ind + terminal2.ind) + nvar + 3 + 2*jointGeneral)
}

#########################################################################
#########################################################################
# (5) Define GH nodes Weights

gh <- statmod::gauss.quad(GHpoints, kind="hermite")
ghNodes = gh$nodes
ghWeights = gh$weights * exp(gh$nodes^2)

############################################################
############################################################
# (6) Compute Gap Times (If Applicable)

if(gapTimes){
	tt11 <- tt11 - tt10
	tt10 <- 0 * tt10
	tt1meta0 <- tt1meta0 - tt0meta0
	tt0meta0 <- 0 * tt0meta0
}

#########################################################################
#########################################################################
# (7) Fill Parameter Vector with User-Defined Values OR Initialize Models

# Check if user entered values for hazard, input 1s if not
if(is.null(init.hazard)) init.hazard <- rep(1, np - nvar - 3 - 2*jointGeneral)
# Check if user entered values for coefficients, input 0s if not
if(is.null(init.B)) init.B <- rep(0, nvar)

# Check lengths of inputs
if(typeof == "Weibull" & length(init.hazard != 2 * (2 + event2.ind + terminal2.ind))){
	stop("init.hazard must have length 6 for weibull for three
	     event types, or length 8 for four event types.")
}else if(typeof == "Splines" & length(init.hazard != (n.knots + 2) * (2 + event2.ind + terminal2.ind))){
		stop("init.hazard must have length (n.knots + 2) * number of event types (3 or 4) for splines.")
}
if(jointGeneral & length(init.Sigma)!=3){
	stop("init.Sigma must have length 3 when jointGeneral = T.\n
	     Order should be: c(frailtySDTerminal1,  frailtySDTerminal2, frailtyCorrelation)")
}else if(!jointGeneral & length(init.Sigma)!=1){
	stop("init.Sigma must have length 1 when jointGeneral = F.")
}
if(length(init.B) != nvar){
	stop("init.B must be the same length as the number of coefficients.")
}

# If initialization desired, replace values
if(initialize){

	# ignore user-supplied initialization values if initialize == T
	init.hazard <- init.hazard*0 + 1

	init.B <- init.B*0

	# recreate time variable in original data set in case of gap times, create new formula
	if(gapTimes){
		initialization.formula <-
			paste("Surv(gapTimes, ", EVENT1, ")",
			      paste(gsub("Surv(.*)","", as.character(formula)), collapse = ""),
			      collapse = "")

		data$gapTimes <- tt11
	}else{
		initialization.formula <- formula
	}
	initialization.formula <- as.formula(initialization.formula)

	# create formula for initialization model 1 (includes terminal 1)
	initialization.formula <- terms(initialization.formula, specials = specials)
	initialization.formula1 <- drop.terms(terms(initialization.formula),
				  survival::untangle.specials(terms(initialization.formula, c("terminal2")), "terminal2", 1:10)$terms,
				  keep.response = T)
	initialization.formula1 <- formula(initialization.formula1)

	# create formula for initialization model 2 (includes terminal 2)
	initialization.formula2 <- drop.terms(terms(initialization.formula),
				  survival::untangle.specials(terms(initialization.formula, c("terminal")), "terminal", 1:10)$terms,
				  keep.response = T)
	initialization.formula2 <- formula(initialization.formula2)
	initialization.formula2 <- sub("terminal2\\(","terminal\\(",initialization.formula2)
	initialization.formula2 <- formula(paste0(initialization.formula2[2:3], collapse = "~"))


	if(!(hazard %in% c("Splines","Splines-per"))){
		mod.joint1<-frailtyPenal(formula = initialization.formula1,
				formula.terminalEvent = formula.terminalEvent,
				jointGeneral = F,
				data = data,
				recurrentAG = !gapTimes,
				hazard = "Weibull",
				RandDist = "LogN",
				maxit = 100, print.times = F)

		# fit initialization model 2
		mod.joint2<-frailtyPenal(formula = initialization.formula2,
				formula.terminalEvent = formula.terminalEvent2,
				jointGeneral = F,
				data = data,
				recurrentAG = !gapTimes,
				hazard = "Weibull",
				RandDist = "LogN",
				maxit = 100, print.times = F)
	}else{
			# fit initialization model 1
			mod.joint1<-frailtyPenal(formula = initialization.formula1,
					formula.terminalEvent = formula.terminalEvent,
					jointGeneral = F,
					data = data,
					kappa=kappa[1:2],
					recurrentAG = !gapTimes,
					hazard = "Splines", RandDist = "LogN",
					n.knots=n.knots,
					maxit = 100, print.times = F)

			# fit initialization model 2
			mod.joint2<-frailtyPenal(formula = initialization.formula2,
					formula.terminalEvent = formula.terminalEvent2,
					jointGeneral = F,
					data = data,
					recurrentAG = !gapTimes,
					n.knots=n.knots,
					kappa=kappa[c(1,3)],
					hazard = "Splines", RandDist = "LogN",
					maxit = 100, print.times = F)

	}
	# Grab initialized values
		# Note: Joint model optimizes on the square root scale
		# for hazard parameters and frailty
		# variance, so we have to square to get to the original scale.
	# Recurrent Hazard

	if(!(hazard %in% c("Splines","Splines-per"))){ #weibull
		init.hazard[1:2]=(mod.joint1$b[1:2]^2+mod.joint2$b[1:2]^2)/2
		init.hazard[3:4]=mod.joint1$b[3:4]^2
		init.hazard[5:6]=mod.joint2$b[3:4]^2
	}else{ #splines 
		init.hazard[1:(2+n.knots)] <- (mod.joint1$b[1:(2+n.knots)]^2 + mod.joint2$b[1:(2+n.knots)]^2)/2
		# average estimates from the two models
		# Terminal 1 Hazard
		init.hazard[(3+n.knots):(4+n.knots*2)] <- mod.joint1$b[(3+n.knots):(4+n.knots*2)]^2
		# Terminal 2 Hazard
		init.hazard[(5+n.knots*2):(6+n.knots*3)] <- mod.joint2$b[(3+n.knots):(4+n.knots*2)]^2
	}
	# Random Effect Variance
	if(!jointGeneral){
		init.Sigma <- (abs(mod.joint1$b[5+n.knots*2]) + abs(mod.joint2$b[5+n.knots*2]))/2
			# average estimates from the two models
	}else{
		init.Sigma[1] <- abs(mod.joint1$b[5+n.knots*2])
		init.Sigma[2] <- abs(mod.joint2$b[5+n.knots*2])
		init.Sigma[3] <- 0 # rho, covariance
	}
	# Alpha
	init.Alpha1 <- mod.joint1$b[6+n.knots*2]
	init.Alpha2 <- mod.joint2$b[6+n.knots*2]

	# Coefficients
	if(noVarEvent[1] == 0){
		# average two estimates
		init.B[1:nbvar[1]] <- (mod.joint1$b[(7+n.knots*2):(6+n.knots*2+nbvar[1])] + mod.joint2$b[(7+n.knots*2):(6+n.knots*2+nbvar[1])])/2
	}
	if(noVarEvent[2] == 0){
		init.B[(1+nbvar[1]):(nbvar[1]+nbvar[2])] <- mod.joint1$b[(7+n.knots*2+nbvar[1]):(6+n.knots*2+nbvar[1]+nbvar[2])]
	}
	if(noVarEvent[4] == 0){
		init.B[(1+nbvar[1]+nbvar[2]):(nbvar[1]+nbvar[2]+nbvar[4])] <- mod.joint2$b[(7+n.knots*2+nbvar[1]):(6+n.knots*2+nbvar[1]+nbvar[4])]
	}
}

# Fill parameter vector
if(!jointGeneral){
	b <- c(sqrt(init.hazard),
	       log(init.Sigma),
	       init.Alpha1, init.Alpha2,
	       init.B)
}else{
	b <- c(sqrt(init.hazard),
	       log(init.Sigma[1:2]), # variance
	       log((init.Sigma[3]+1)/(1-init.Sigma[3])), # rho transformed using scale-logit
	       init.Alpha1, init.Alpha2,
	       init.B)
}
# save a copy of the starting value for the output
start.b <- b

if(length(b)!=np) stop("Parameter vector not the correct length.")

if(crossVal==T & hazard %in% c('Splines','Splines-per')){
	if(event2.ind==0){
	xx=as.character(initialization.formula1)
	newf=paste(xx[[2]],'~1',sep='')

	mm1=frailtyPenal(as.formula(newf), data=data,cross.validation = T,
				n.knots = n.knots, kappa = 10,print.times = F)
	mm2=frailtyPenal(Surv(tt1dc,terminal1)~1,cross.validation = T,
				n.knots = n.knots, kappa = 10,
				data=data.frame(tt1dc=tt1dc,terminal1=terminal1),
				print.times = F)
	mm3=frailtyPenal(Surv(tt1dc,terminal2)~1,cross.validation = T,
				n.knots = n.knots, kappa = 10,
				data=data.frame(tt1dc=tt1dc,terminal2=terminal2),
				print.times = F)
	kappa=c(mm1$kappa,mm2$kappa,0,mm3$kappa)
	}
	if(initialize){
		betah1 = mm1$b[1:(2+n.knots)]^2
		betah2 = mm2$b[1:(2+n.knots)]^2
		betah3 = mm3$b[1:(2+n.knots)]^2
	b <- c(sqrt(c(betah1,betah2,betah3)),
	       log(init.Sigma),
	       init.Alpha1, init.Alpha2,
	       init.B)
	}
}

############################################################
############################################################
# (8) Check Dimensions of all variables for debugging
controls = c(maxit = maxit[1], # [1]
	 initialize = initialize, # [2]
	 typeof = typeof, # [3]
	 equidistant=equidistant, #[4]
	 irep = !crossVal, # [5] irep
	 gapTimes = gapTimes, # [6] ag0
	 nbIntervEvent = 0, # [7] nbIntervEvent
	 n.knots = n.knots, # [8]
	 event2.ind = event2.ind, # [9]
	 terminal2.ind = terminal2.ind, # [10]
	 GHpoints = GHpoints, # [11]
	 jointGeneral = as.integer(jointGeneral)) # [12] typeJoint0
if(length(controls) != 12) stop("Length of 'controls' not 12.")


nobsEvent = c(length(event1),
	  length(terminal1),
	  length(event2))
if(length(nobsEvent)!= 3) stop("Length of 'nobsEvent' not 3.")


if(length(kappa) != 4) stop("Length of 'kappa' not 4.")


### Recurrent 1
if(length(tt10) != nobsEvent[1] | length(tt10) !=  nobsEvent[1] | length(event1) !=  nobsEvent[1] | length(group1) !=  nobsEvent[1]){
	stop("Length of tt00, tt10, event1, group1 not nobsEvent[1]")
}


### Recurrent 2
if(length(tt0meta0) !=  nobsEvent[3]  | length(tt1meta0) !=  nobsEvent[3]  | length(group2) !=  nobsEvent[3]  | length(event2) !=  nobsEvent[3] ){
	stop("Length of tt0meta0, tt1meta0, group2, event2 not nobsEvent[3]")
}


if(length(event1) != nobsEvent[1]) stop("length(event1) != nobsEvent[1]")
if(length(event2) != nobsEvent[3]) stop("length(event2) != nobsEvent[3]")
if(length(group1) != nobsEvent[1]) stop("length(group1) != nobsEvent[1]")
if(length(group2) != nobsEvent[3]) stop("length(group2) != nobsEvent[3]")
if(length(groupdc) != nobsEvent[2]) stop("length(groupdc) != nobsEvent[2]")
if(max(group1) != nobsEvent[2]) stop("max(group1) != nobsEvent[2]")
if(max(groupdc) != nobsEvent[2]) stop("max(groupdc) != nobsEvent[2]")
if(length(tt1dc) != nobsEvent[2]) stop("length(tt1dc) != nobsEvent[2]")

### Terminal Events
if(length(terminal1) != nobsEvent[2]){
	stop("length(terminal1) != nobsEvent[2]")
}
if(length(terminal2) != nobsEvent[2]){
	stop("length(terminal2) != nobsEvent[2]")
}

if(length(nbvar) != 4) stop("length(nbvar) != 4")


if(all(dim(modelmatrix1) != c(nobsEvent[1],nbvar[1]))) stop("all(dim(modelmatrix1) != c(nobsEvent[1],nbvar[1]))")
if(all(dim(modelmatrix2) != c(nobsEvent[2],nbvar[2]))) stop("all(dim(modelmatrix2) != c(nobsEvent[3],nbvar[2]))")
if(all(dim(modelmatrix3) != c(nobsEvent[3],nbvar[3]))) stop("all(dim(modelmatrix3) != c(nobsEvent[2],nbvar[3]))")
if(all(dim(modelmatrix4) != c(nobsEvent[2],nbvar[4]))) stop("all(dim(modelmatrix4) != c(nobsEvent[3],nbvar[4]))")


if(length(noVarEvent) != 4) stop("length(noVarEvent) != 4")


if(any(is.na(modelmatrix1))|any(is.na(modelmatrix2))|any(is.na(modelmatrix3))|any(is.na(modelmatrix4))){
	stop("NA values among covariates. Reconfigure Data.")
}
######################################################################################################
# (9) Send to Fortran for optimization 

    ans <- .Fortran(C_joint_competing,
                controls = as.integer(controls),
                nobsEvent = as.integer(nobsEvent), #nobsEvent
                k0 = as.double(kappa),

                # Data Arguments
                tt00 = as.double(tt10),
                tt10 = as.double(tt11),
                tt0meta0 = as.double(tt0meta0),
                tt1meta0 = as.double(tt1meta0),
                ic0 = as.integer(event1), #ic0
                icmeta0 = as.integer(event2), #icmeta0

                groupe0 = as.integer(group1),#groupe0
                groupe0meta = as.integer(group2),#groupe0meta
                groupe0dc = as.integer(groupdc),#groupe0dc

                tt0dc0 = as.double(0*tt1dc), #tt0dc0
                tt1dc0 = as.double(tt1dc), #tt1dc0
                icdc0 = as.integer(terminal1), #icdc0
                icdc20 = as.integer(terminal2), #icdc20

                nbvar = as.integer(nbvar), #nbvar

                vax0 = as.double(modelmatrix1),
                vaxdc0 = as.double(modelmatrix2), #vaxdc0
                vaxmeta0 = as.double(modelmatrix3),
                vaxdc20 = as.double(modelmatrix4), #vaxdc20

                noVarEvent = as.integer(noVarEvent), #noVarEvent

                # Parameter Information
                np=as.integer(np),
                b=as.double(b),
                H_hessOut=as.double(matrix(0,nrow=np,ncol=np)),
                HIHOut=as.double(matrix(0,nrow=np,ncol=np)),
                resOut=as.double(0),
                LCV=as.double(rep(0,2)),
                critCV=as.integer(rep(0,3)),
                mtEvent = as.integer(rep(100,4)), #mtEvent
                mt1Event = as.integer(rep(100,4)), #mt1Event

                # Survival and Hazard Function Fits
                x1=as.double(rep(0,100)),                  #x1Out
                lam=as.double(matrix(0,nrow=100,ncol=3)),  #lamOut
                xSu1=as.double(rep(0,100)),                #xSu1
                surv=as.double(matrix(0,nrow=100,ncol=3)), #suOut
                x2=as.double(rep(0,100)),                  #x2Out
                lam2=as.double(matrix(0,nrow=100,ncol=3)), #lam2Out
                xSu2=as.double(rep(0,100)),                #xSu2
                surv2=as.double(matrix(0,nrow=100,ncol=3)),#su2Out
                x3=as.double(rep(0,100)),                  #x3Out
                lam3=as.double(matrix(0,nrow=100,ncol=3)), #lam3out
                xSu3=as.double(rep(0,100)),                #xSu3
                surv3=as.double(matrix(0,nrow=100,ncol=3)),#su3Out
                x4=as.double(rep(0,100)),                  #x4Out
				lam4=as.double(matrix(0,nrow=100,ncol=3)), #lam4Out
				xSu4=as.double(rep(0,100)),                #xSu4
				surv4=as.double(matrix(0,nrow=100,ncol=3)),#su4Out
                ni=as.integer(0),
				constraintfrailty=as.integer(constraint.frailty),
                cptEvent=as.integer(rep(0,4)),
                ResMartingaleEvent=as.double(matrix(0,nrow=nobsEvent[2],ncol=3)),
                frailtyEstimates=as.double(matrix(0,nrow=nobsEvent[2],ncol=5)),
                linearpred=as.double(rep(0,nobsEvent[1])),
                linearpreddc=as.double(rep(0,nobsEvent[2])),
                linearpredM=as.double(rep(0,nobsEvent[3])),
                linearpreddc2=as.double(rep(0,nobsEvent[2])),
                ziOut1=as.double(rep(0,controls[8]+6)),
                ziOutdc=as.double(rep(0,controls[8]+6)),
                ziOutmeta=as.double(rep(0,controls[8]+6)),
                time=as.double(rep(0,controls[7]+1)),
                timedc=as.double(rep(0,controls[7]+1)),
                timeM=as.double(rep(0,controls[7]+1)),
                ghNodes = as.double(ghNodes),
                ghWeights = as.double(ghWeights),
                tolerance0 = as.double(tolerance)
    )
 if(ans$critCV[2] != 1){
	cat('Model did not converge.','\n')
	ans<-NULL
	return(ans)
 }
 ans$call=call
######################################################################################################
######################################################################################################
# (10) Format Model Summary Tables
 if(jointGeneral == F & hazard == "Weibull"){
 	f <- function(b){
 		c(b[1:6]^2,#exp(b[1:6])^2,
 		  ifelse(constraint.frailty,exp(b[7]),abs(b[7])),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       ifelse(constraint.frailty,exp(b[7]),2*b[7]),
 		       rep(1,nvar + 2)))
 	}
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Sigma",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
	nhaz = 6 
	ans$varH.Raw <- matrix(ans$H_hessOut, nrow = np, ncol = np)
	ans$varH.Estimate <- f.prime(ans$b) %*% ans$varH.Raw %*% f.prime(ans$b)
	Raw.SE<-NULL
	Raw<-NULL
	Estimate<-NULL
	Estimate.SE<-NULL
	ans$summary.table <- data.frame(
		Parameter = Parameter,
		Raw = ans$b,
		Raw.SE = sqrt(diag(ans$varH.Raw)),
		Estimate = f(ans$b),
		Estimate.SE = sqrt(diag(ans$varH.Estimate)),
		LB95 = f(ans$b - 2*sqrt(diag(ans$varH.Raw))),
		UB95 = f(ans$b  + 2*sqrt(diag(ans$varH.Raw))),
		p = 2*pnorm(q = -abs(ans$b), mean = 0,
					sd = sqrt(diag(ans$varH.Raw))),
		H0 = paste(Parameter, " = ", f(rep(0, np)))
	)
 }else if(jointGeneral == T & hazard == "Weibull"){
 	f <- function(b){
 		c(b[1:6]^2,#exp(b[1:6])^2,
 		  exp(b[7:8]),
 		  (exp(b[9]) - 1)/(exp(b[9]) + 1),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       exp(b[7:8]),
 		       2*exp(2*b[9])/(1+exp(b[9]))^2,
 		       rep(1,nvar+2)))
 	}
	nhaz=6
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Terminal1: Sigma", "Terminal2: Sigma", "Rho",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
	ans$varH.Raw <- matrix(ans$HIHOut, nrow = np, ncol = np)
	ans$varH.Estimate <- f.prime(ans$b) %*% ans$varH.Raw %*% f.prime(ans$b)
	Raw.SE<-NULL
	Raw<-NULL
	Estimate<-NULL
	Estimate.SE<-NULL
	ans$summary.table <- data.frame(
		Parameter = Parameter,
		Raw = ans$b,
		Raw.SE = sqrt(diag(ans$varH.Raw)),
		Estimate = f(ans$b),
		Estimate.SE = sqrt(diag(ans$varH.Estimate)),
		LB95 = f(ans$b - 2*sqrt(diag(ans$varH.Raw))),
		UB95 = f(ans$b  + 2*sqrt(diag(ans$varH.Raw))),
		p = 2*pnorm(q = -abs(ans$b), mean = 0, 
					sd = sqrt(diag(ans$varH.Raw))),
		H0 = paste(Parameter, " = ", f(rep(0, np)))
	)
 }else if(jointGeneral == F & hazard %in% c("Splines","Splines-per")){
	nhaz=(2+n.knots)*(2+event2.ind+terminal2.ind)
	f <- function(b){c(exp(b[nhaz+1]),
 		    b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(exp(b[nhaz+1]),
 		       rep(1,nvar + 2)))
	}	
 	Parameter = c(
 		  "Sigma",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
	#ans$varH.Raw <- matrix(ans$HIHOut, nrow = np, ncol = np)
	ans$varH.Raw <- matrix(ans$H_hessOut, nrow = np, ncol = np)
	ans$varH.Estimate <- f.prime(ans$b) %*% ans$varH.Raw[(nhaz+1):np,(nhaz+1):np] %*% f.prime(ans$b)
	Raw.SE<-NULL
	Raw<-NULL
	Estimate<-NULL
	Estimate.SE<-NULL
	ans$summary.table <- data.frame(
		Parameter = Parameter,
		Raw = ans$b[(nhaz+1):np],
		Raw.SE = sqrt(diag(ans$varH.Raw))[(nhaz+1):np],
		Estimate = f(ans$b),
		Estimate.SE = sqrt(diag(ans$varH.Estimate)),
		LB95 = f(ans$b - 2*sqrt(diag(ans$varH.Raw))),
		UB95 = f(ans$b  + 2*sqrt(diag(ans$varH.Raw))),
		p = 2*pnorm(q = -abs(ans$b[(nhaz+1):np]), mean = 0,
		sd = sqrt(diag(ans$varH.Raw))[(nhaz+1):np]),
		H0 = paste(Parameter, " = ", f(rep(0, np)))
	)
 }else if(jointGeneral == T & hazard %in% c("Splines","Splines-per")){
	nhaz=(2+n.knots)*(2+event2.ind+terminal2.ind)
	f <- function(b){c(ifelse(constraint.frailty,exp(b[nhaz+1]),abs(b[nhaz+1])),
 		    b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(ifelse(constraint.frailty,exp(b[nhaz+1]),2*b[nhaz+1]),
 		       rep(1,nvar + 2)))
	}	
 	Parameter = c(
 		  "Terminal1: Sigma", "Terminal2: Sigma", "Rho",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
	ans$varH.Raw <- matrix(ans$HIHOut, nrow = np, ncol = np)
	ans$varH.Estimate <- f.prime(ans$b) %*% ans$varH.Raw[(nhaz+1):np,(nhaz+1):np] %*% f.prime(ans$b)
	Raw.SE<-NULL
	Raw<-NULL
	Estimate<-NULL
	Estimate.SE<-NULL
	ans$summary.table <- data.frame(
		Parameter = Parameter,
		Raw = ans$b[(nhaz+1):np],
		Raw.SE = sqrt(diag(ans$varH.Raw))[(nhaz+1):np],
		Estimate = f(ans$b),
		Estimate.SE = sqrt(diag(ans$varH.Estimate)),
		LB95 = f(ans$b - 2*sqrt(diag(ans$varH.Raw))),
		UB95 = f(ans$b  + 2*sqrt(diag(ans$varH.Raw))),
		p = 2*pnorm(q = -abs(ans$b[(nhaz+1):np]), mean = 0, 
					sd = sqrt(diag(ans$varH.Raw))[(nhaz+1):np]),
		H0 = paste(Parameter, " = ", f(rep(0, np)))
	)
 }

######################################################################################################
######################################################################################################
# (11) Format Initialization Tables
 ans$initialization$b <- start.b
 if(initialize & FALSE){
 	ans$initialization$joint1 <- mod.joint1
 	ans$initialization$joint2 <- mod.joint2

	# Extract Transformed Parameters
 	f1 <- function(b, i=3){
 		c(b[1:(length(b)-nvar+nbvar[i]-2)]^2,
 		  b[(length(b)-nvar+nbvar[i]-1):length(b)])
 	}
 	f1.prime <- function(b, i = 3){
 		diag(c(2*b[1:(length(b)-nvar+nbvar[i]-2)],
 		       rep(1,nvar-nbvar[i]+2)))
 	}
 	Parameters1 = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Sigma", "Terminal1: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)))

 	Parameters2 = c("Recurrent: Shape", "Recurrent: Scale",
 		    "Terminal2: Shape", "Terminal2: Scale",
 		    "Sigma", "Terminal2: Alpha",
 		    paste0("Recurrent: ",colnames(modelmatrix1)),
 		    paste0("Terminal2: ",colnames(modelmatrix4)))

 	ans$initialization$varH.Estimate1 <-
 		f1.prime(ans$initialization$joint1$b) %*%
 		ans$initialization$joint1$varHtotal %*%
 		f1.prime(ans$initialization$joint1$b)
 	ans$initialization$varH.Estimate2 <-
 		f1.prime(ans$initialization$joint2$b, i=4) %*%
 		ans$initialization$joint2$varHtotal %*%
 		f1.prime(ans$initialization$joint2$b, i=4)
 	ans$initialization$summary.table1 <- data.frame(
 		Parameter = Parameters1,
 		Raw = ans$initialization$joint1$b,
 		Raw.SE = sqrt(diag(ans$initialization$joint1$varHtotal)),
 		Estimate = f1(ans$initialization$joint1$b),
 		Estimate.SE = sqrt(diag(ans$initialization$varH.Estimate1)),
 		p = 2*pnorm(q = -abs(ans$initialization$joint1$b), mean = 0,
					sd =sqrt(diag(ans$initialization$joint1$varHtotal)))
 	)
 	ans$initialization$summary.table2 <- data.frame(
 		Parameter = Parameters2,
 		Raw = ans$initialization$joint2$b,
 		Raw.SE = sqrt(diag(ans$initialization$joint2$varHtotal)),
 		Estimate = f1(ans$initialization$joint2$b, i=4),
 		Estimate.SE = sqrt(diag(ans$initialization$varH.Estimate2)),
 		LB95 = f1(ans$initialization$joint2$b, i=4) - 2*sqrt(diag(ans$initialization$varH.Estimate2)),
 		UB95 = f1(ans$initialization$joint2$b, i=4) + 2*sqrt(diag(ans$initialization$varH.Estimate2)),
 		p = 2*pnorm(q = -abs(ans$initialization$joint2$b), mean = 0,
					sd =sqrt(diag(ans$initialization$joint2$varHtotal)))
 	)
 }

 if(typeof==0){
	ans$logLikPenal=ans$resOut
	ans$LCV=ans$LCV[1]
	ans$resOut<-NULL 
 }else{
	ans$logLik=ans$resOut
	ans$AIC=ans$LCV[2]
	ans$LCV<-NULL
	ans$resOut<-NULL 
	ans$ziOut1<-NULL 
	ans$ziOutdc<-NULL
	ans$k0<-NULL
 }




 ans$lam=matrix(ans$lam,100,3)
 ans$lam2=matrix(ans$lam2,100,3)
 ans$lam3=matrix(ans$lam3,100,3)

 ans$surv=matrix(ans$surv,100,3)
 ans$surv2=matrix(ans$surv2,100,3)
 ans$surv3=matrix(ans$surv3,100,3)


 ans$n = length(unique(ans$groupe0))
 ans$nevts = c(sum(ans$ic0),sum(ans$icdc0),sum(ans$icdc20))

 names(ans)[which(names(ans)=="tolerance0")]="tolerance"
 ans$ziOutmeta<-NULL
 ans$time<-NULL
 ans$timedc<-NULL
 ans$timeM<-NULL
 ans$cptEvent<-NULL
 ans$ResMartingaleEvent<-NULL
 ans$frailtyEstimates<-NULL
 ans$linearpred<-NULL
 ans$linearpreddc<-NULL
 ans$linearpredM<-NULL
 ans$linearpreddc2<-NULL
 ans$x4<-NULL
 ans$surv4<-NULL
 ans$lam4<-NULL
 ans$xSu4<-NULL
 ans$tt00<-NULL
 ans$tt10<-NULL
 ans$tt0meta0<-NULL
 ans$tt1meta0<-NULL
 ans$ic0<-NULL
 ans$icmeta0<-NULL
 ans$groupe0<-NULL
 ans$groupe0meta<-NULL
 ans$groupe0dc<-NULL
 ans$tt0dc0<-NULL
 ans$tt1dc0<-NULL
 ans$icdc0<-NULL
 ans$icdc20<-NULL
 ans$nbvar<-NULL
 ans$vax0<-NULL
 ans$vaxdc0<-NULL
 ans$vaxmeta0<-NULL
 ans$vaxdc20<-NULL
 ans$nobsEvent<-NULL
 ans$mtEvent<-NULL
 ans$mt1Event<-NULL
 ans$varH.Raw<-NULL
 ans$varH.Estimate<-NULL
 ans$initialization<-NULL
 class(ans)<-"jointRecCompet"
return(ans)
}