##' Short summary of the random effects parameters, the fixed treatment
##' effects, and the surrogacy evaluation criteria estimated from a joint surrogate mediation model
##'
##' This function returns the estimate of the coefficients of the model, their standard error and the
##' associated p-values of the Wald test for the joint surrogate model, also hazard ratios (HR) and their
##' confidence intervals for the fixed treatment effects. It also displays summary of the surrogacy measure \eqn{R(t)}
##' and of the natural direct, indirect and total effect.
##'
##' @aliases summary.jointSurroMed print.summary.jointSurroMed
##' @usage \method{summary}{jointSurroMed}(object,d=4,len=3,n=3,...)
##'
##' @param object An object inheriting from \code{jointSurroMed} class.
##' @param d The desired number of digits after the decimal point for parameters.
##' The maximum of 4 digits is required for the estimates. Default of 3 digits is used.
##' @param len The desired number of digits after the decimal point for p-value and convergence
##' criteria. Default of 4 digits is used.
##' @param n The number of time points to be used in the results of the differents function
##' related to the mediation analysis: \eqn{g(s)}, \eqn{R(t)} and the direct, indirect and total
##' effect. The provided value should be between 1 and 20. Default is 3.
##' @param ... other unused arguments.
##'
##' @return For the variances parameters of the random effects, it prints the estimate of
##' the coefficients with their standard error, Z-statistics and p-values
##' of the Wald test. For the fixed treatment effects, it also prints HR and its confidence
##' intervals for each covariate.
##' For the surrogacy assessment, prints \code{n} value of the estimation function \eqn{g(s)} and \eqn{R(t)}.
##' Also prints the values of the estimated direct, indirect and total effects.

##' The remaining displayed information concern the convergence characteristics and
##' include the penalized marginal log-likelihood, the number of iterations, the LCV and the convergence criteria.
##' @seealso \code{\link{jointSurroPenal}}
##' @export
##' @keywords methods

"summary.jointSurroMed"<-
function(object, d = 4, len = 3,n=3,...){
  x <- object
  int.method.kt = 0
  if (!inherits(x, "jointSurroMed"))
    stop("Object must be of class 'jointSurroMed'")

  if(n>10){
    stop("The argument 'n' should be less than 10.")
  }


  cat("Results of a joint surrogate mediation model using a penalized likelihood\non the baseline hazard functions.\n")
  cat(" ", "\n")

  coef <- data.frame(x$Coefficients)

  beta <- coef[(1 : (nrow(coef)-2)),1:2]
  names(beta)[2] <- "SE"
  beta$"z" <- round(beta$Estimate/beta$SE,len)
  beta$P <- signif(1 - pchisq((beta$Estimate/beta$SE)^2, 1), 5)
  beta$" " <- ifelse(beta$P < 0.001,"***",ifelse(beta$P < 0.01,"**",
                                                 ifelse(beta$P < 0.05,"*",ifelse(beta$P < 0.1,"."," "))))
  beta$P <- formatC(beta$P, d, format = "g")
  beta$Estimate <- round(beta$Estimate,min(4,len))
  names(beta)[2] <-"Std Error"

  # travail des P
  p <- NULL
  p <- ifelse(as.numeric(beta$P) < 10^-10, "< e-10", beta$P)
  beta$P <- p

  cat("Estimates for variances parameters of the random effects", "\n")
  rownames(beta)[(nrow(beta) - 4)] <- "sigma2_S"
  rownames(beta)[(nrow(beta) - 3)] <- "sigma2_T"
  rownames(beta)[(nrow(beta) - 2)] <- "sigma_ST"
  print(beta[1:(nrow(beta) - 2),])

  beta2 <- beta[((nrow(beta) - 1) : nrow(beta)),]
  beta2$Estimate <- round(beta2$Estimate,min(4,len))

  cat(" ", "\n")
  cat("Estimates for the fixed treatment effects", "\n")
  print(beta2)

  cat("---","\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

  cat(" ", "\n")
  cat(" ","\n")

  cat("Hazard ratios (HR) and confidence intervals for the fixed treatment effects", "\n")
  HR <- round(exp(coef[((nrow(coef) - 3) : (nrow(coef)-2)),-2]), len)
  names(HR)[1] <- c("exp(coef)")
  print(HR)

  coef <- rbind(coef,coef[nrow(coef)-1,])
  coef[nrow(coef),c(3,4)] <- object$R2.boot[-1]
  coef[nrow(coef),1] <- object$R2.boot[1]
  coef[nrow(coef),2] <- NA
  rownames(coef)[nrow(coef)] <- "R2.boot"


  validation <- coef[c(nrow(coef) - 1, nrow(coef) - 2,nrow(coef)),]
  validation[,2] <- as.character(validation[,2])
  if(x$type.joint == 1) {
    validation[1,2] <- "--"
  }else{
    validation[1,2] <- as.character(round(as.numeric(validation[1,2]), len))
  }
  validation[3,2] <- "--"
  validation[,1] <- round(validation[,1], len)
  validation[,3] <- round(validation[,3], len)
  validation[,4] <- round(validation[,4], len)
  validation[2,2] <- as.character(round(as.numeric(validation[2,2]), len))
  if(int.method.kt == 1){
    validation[1,1] <- jointSurroTKendall(object = object, int.method = 1)
    # validation[1,1] <- jointSurroTKendall(theta = object$Coefficients["theta",1],
    #                                       gamma = ifelse(is.na(object$Coefficients["gamma",1]), 0, object$Coefficients["gamma",1]),
    #                                       alpha = ifelse(is.na(object$Coefficients["alpha",1]), 1, object$Coefficients["alpha",1]),
    #                                       zeta = ifelse(is.na(object$Coefficients["zeta",1]), 1, object$Coefficients["zeta",1]),
    #                                       nb.gh = nb.gh, ui = ifelse(is.na(object$Coefficients["gamma",1]), 0, 1), int.method = 1)
  }

  names(validation)[2] <- "Std Error"

  # More precision on the surrogacy evaluation
  validation2 <- data.frame(matrix(rep(NA,18), nrow = 3, ncol = 6))
  names(validation2) <- c("Level", names(validation), "Strength")
  rownames(validation2) <- rownames(validation)
  validation2[,1] <- c("Individual", "Trial", "Trial")
  validation2[,2:5] <- validation
  validation2[2,6] <- ifelse(validation2[2,4] <= 0.49,"Low",ifelse(validation2[2,4]<0.72,
                                                                   "Medium","High"))
  validation2[3,6] <- ifelse(validation2[3,4] <= 0.49,"Low",ifelse(validation2[2,4]<0.72,
                                                                   "Medium","High"))
  validation2[1,6] <- " "

  cat(" ", "\n")
  cat("Individual and trial level associations", "\n")
  print(validation2)
  cat("---","\n")
  cat("Correlation strength: <= 0.49 'Low'; ]0.49 - 0.72[ 'Medium'; >= 0.72 'High' ","\n")
  cat("---","\n")


  ## mediation ...
  mediation<-object$mediation
  n<-min(10,n,length(mediation$data.rt$Time),length(mediation$data.g$s))
  ## Gamma
  times<-mediation$data.g$s
  g<-mediation$data.g$g
  upp<-mediation$data.g$upper
  low<-mediation$data.g$lower


  pos<-(1:n)*floor((length(times)/n))

  times<-round(times[pos],min(4,len))
  g<-round(g[pos],min(4,len))
  upp<-round(upp[pos],min(4,len))
  low<-round(low[pos],min(4,len))

  Ci = sapply(seq_along(low),function(i){
    paste("[",low[i],";",upp[i],"]",sep="")
  })

  cat(" ", "\n")
  cat("Estimated function g at", n, "times of occurence of the surrogate\n")
  print(data.frame("Time"=times,"g"=g,"CI 95"=Ci))

  ## R(t) and NEFF ......
  times<-round(mediation$data.rt$Time,min(4,len))
  rt<-round(mediation$data.rt$Rt,min(4,len))
  nie<-round(mediation$data.rt$NIE,min(4,len))
  nde<-round(mediation$data.rt$NDE,min(4,len))
  te<-round(mediation$data.rt$TE,min(4,len))

  pos<-(1:n)*floor((length(times)/n))

  times<-times[pos]
  rt<-rt[pos]
  nie<-nie[pos]
  nde<-nde[pos]
  te<-te[pos]
  cat(" ", "\n")
  cat("Estimated function R(t), natural direct, indirect and total effect at",
      n, "time points \n")

  if(length(mediation)==5){
    #no confidence bands available
    print(data.frame("Time"=times,"Rt"=rt,"Total"=te,
                     "Direct"=nde,"Indirect"=nie))
  }
  if(length(mediation)==9){
    nie.upp<-round(mediation$NIE.ci$upper[pos],min(4,len))
    nie.low<-round(mediation$NIE.ci$lower[pos],min(4,len))
    nde.upp<-round(mediation$NDE.ci$upper[pos],min(4,len))
    nde.low<-round(mediation$NDE.ci$lower[pos],min(4,len))
    te.upp<-round(mediation$TE.ci$upper[pos],min(4,len))
    te.low<-round(mediation$TE.ci$lower[pos],min(4,len))
    rt.upp<-round(mediation$Rt.ci$upper[pos],min(4,len))
    rt.low<-round(mediation$Rt.ci$low[pos],min(4,len))

    nie.ci = sapply(seq_along(nie.upp),function(i){
      paste("[",nie.low[i],";",nie.upp[i],"]",sep="")
    })
    nde.ci = sapply(seq_along(nde.upp),function(i){
      paste("[",nde.low[i],";",nde.upp[i],"]",sep="")
    })
    te.ci = sapply(seq_along(te.upp),function(i){
      paste("[",te.low[i],";",te.upp[i],"]",sep="")
    })
    rt.ci = sapply(seq_along(rt.upp),function(i){
      paste("[",rt.low[i],";",rt.upp[i],"]",sep="")
    })
    print(data.frame("Time"=times,"Rt"=rt,"CI 95 Rt"=rt.ci,
                     "TE"=te,"CI 95 TE"=te.ci,
                     "NDE"=nde,"CI 95 NDE"=nde.ci,
                     "NIE"=nie,"CI 95 NIE"=nie.ci))
  }


  cat(" ", "\n")
  cat("Convergence parameters", "\n")
  cat(c("Penalized marginal log-likelihood = ", round(object$loglikPenal, len)), "\n")
  cat(c("Number of iterations = ", object$n.iter),"\n")
  cat("LCV = the approximate likelihood cross-validation criterion", "\n")
  cat(c("      in the semi parametrical case     = ", round(object$LCV, len)),"\n")
  cat("Convergence criteria:", "\n")
  EPS <- formatC(object$EPS, d, format = "g")
  cat(c("  parameters = ",EPS[1], "likelihood = ", EPS[2], "gradient = ", EPS[3]), "\n")
}
