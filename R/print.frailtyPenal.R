
"print.frailtyPenal" <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        if (x$type == "counting") {
            cat("\n      left truncated structure used")
        }
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" frailtyPenal failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- x$coef
    nvar <- length(x$coef)
 
    if (is.null(coef))
      {
        x$varH<-matrix(x$varH) 
        x$varHIH<-matrix(x$varHIH)
      }
    
    if (!is.null(coef)) 
      { 
        seH <- sqrt(diag(x$varH))[-1]
        seHIH <- sqrt(diag(x$varHIH))[-1]
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
        pchisq((coef/seH)^2, 1), digits - 1))
        cat("\n")
        if (!is.null(x$theta)) 
          {
          cat("  Shared Gamma Frailty model parameter estimates ", 
            "\n")
          cat("  using a Penalized Likelihood on the hazard function", 
            "\n")
        if (x$n.strat>1)  
          cat("  (Stratification structure used)", "\n")         
          }
        else {
          cat("  Cox proportional hazards model parameter estimates ", 
            "\n")
          cat("  using a Penalized Likelihood on the hazard function", 
            "\n")
        if (x$n.strat>1)  
          cat("  (Stratification structure used)", "\n")         
          }
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "SE coef (H)", "SE coef (HIH)", "z", "p"))
        cat("\n")
        prmatrix(tmp)
        cat("\n")
      } 
        if (!is.null(x$theta)) {
          tetha <- x$theta
          temp <- diag(x$varH)[1]
          seH <- sqrt(((2 * (tetha^0.5))^2) * temp)
          temp <- diag(x$varHIH)[1]
          seHIH <- sqrt(((2 * (tetha^0.5))^2) * temp)
         cat("    Frailty parameter, Theta:", tetha, "(SE (H):", 
            seH, ")", "(SE (HIH):", seHIH, ")", "\n")
        }
    cat(" \n")
    cat(paste("    penalized marginal log-likelihood =", round(x$logVerComPenal, 
        2)))
    cat("\n")
    cat("    n=", x$n)
    if (length(x$na.action)) 
        cat("  (", length(x$na.action), " observation deleted due to missing) \n")
    else cat("\n")
    cat("    n events=", x$n.event, " n groups=", x$groups)
    cat( "\n")
    cat("    number of iterations: ", x$n.iter)
    cat("\n")
    cat("    Exact number of knots used: ", x$n.knots, "\n")

    if (!x$cross.Val)
      {
       cat("    Value of the smoothing parameter: ", x$kappa[1])
       if (x$n.strat==2)
        cat(" ", x$kappa[2])
      }

    if (x$cross.Val)
      {
       if (is.null(x$theta))
         cat("    Smoothing parameter estimated by Cross validation: ", x$kappa[1])  
       else 
         {
           cat("    Best smoothing parameter estimated by")
           cat("\n")
           cat("       an approximated Cross validation: ", x$kappa[1])
         } 
      }
    
    cat(", DoF: ", formatC(-x$DoF, format="f",dig=2))


    cat("\n")
    invisible()
}
