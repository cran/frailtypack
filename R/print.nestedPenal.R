
"print.nestedPenal" <-
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
        seH <- sqrt(diag(x$varH))[-c(1,2)]
        seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
        pchisq((coef/seH)^2, 1), digits - 1))

        cat("\n")
        cat(" Nested Frailty model parameter estimates using ", 
            "\n")
        cat(" a Penalized Likelihood on the hazard functions", 
            "\n")
        if (x$n.strat>1)  
        cat(" (Stratification structure used)", "\n") 

        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "SE coef (H)", "SE coef (HIH)", "z", "p"))
        cat("\n")
        prmatrix(tmp)
        cat("\n")
      } 

      alpha <- x$alpha
      temp <- diag(x$varH)[1]
      seH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)
      temp <- diag(x$varHIH)[1]
      seHIH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)

      cat("  Frailty parameters: \n")
      cat("   alpha  (group effect): ", alpha, " (SE(H):", 
            seH.alpha, ")", " (SE(HIH):", seHIH.alpha, ")", sep="", "\n")

      eta <- x$eta
      temp <- diag(x$varH)[2]
      seH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)
      temp <- diag(x$varHIH)[2]
      seHIH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)

      cat("   eta (subgroup effect): ", eta, " (SE(H):", 
            seH.eta, ")", " (SE(HIH):", seHIH.eta, ")", sep="", "\n")

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

    cat("    Value of the smoothing parameter: ", x$kappa[1])
    
    cat("\n")
    invisible()
}