#' summary of parameter estimates of a joint frailty model
#'
#' This function returns hazard rations (HR) and its confidence intervals.
#'
#'
#' @aliases summary.jointPenal print.summary.jointPenal
#' @usage \method{summary}{jointPenal}(object, level = 0.95, len = 6, d = 2,
#' lab="HR", ...)
#' @param object output from a call to frailtyPenal for joint models
#' @param level significance level of confidence interval. Default is 95\%.
#' @param d the desired number of digits after the decimal point. Default of 6
#' digits is used.
#' @param len the total field width. Default is 6.
#' @param lab label of printed results.
#' @param \dots other unused arguments.
#' @return Prints HR and its confidence intervals for each covariate.
#' Confidence level is allowed (level argument).
#' @seealso \code{\link{frailtyPenal}}
#' @keywords methods
#' @export
#' @examples
#' \dontrun{
#'
#' data(readmission)
#'
#' #-- gap-time
#' modJoint.gap <- frailtyPenal(
#'   Surv(time, event) ~ cluster(id) + sex + dukes +
#'     charlson + terminal(death),
#'   formula.terminalEvent = ~ sex + dukes + charlson,
#'   data = readmission, n.knots = 14, kappa = c(9.55e+9, 1.41e+12)
#' )
#'
#' #-- calendar time
#' modJoint.calendar <- frailtyPenal(
#'   Surv(t.start, t.stop, event) ~ cluster(id) +
#'     sex + dukes + charlson + terminal(death),
#'   formula.terminalEvent = ~ sex + dukes + charlson,
#'   data = readmission, n.knots = 10, kappa = c(9.55e+9, 1.41e+12), recurrentAG = TRUE
#' )
#'
#' #-- It takes around 1 minute to converge
#'
#' summary(modJoint.gap)
#' summary(modJoint.calendar)
#' }
#'
"summary.jointPenal" <-
  function(object, level = .95, len = 6, d = 2, lab = "HR", ...) {
    if (is.null(object$family)) {
      x <- object
      if (!inherits(x, "jointPenal")) {
        stop("Object must be of class 'frailtyPenal'")
      }

      z <- abs(qnorm((1 - level) / 2))
      co <- x$coef
      indic_alpha <- x$indic_alpha
      if (indic_alpha == 1 || x$joint.clust == 2) {
        se <- sqrt(diag(x$varHIH))[-c(1, 2)]
      } else {
        se <- sqrt(diag(x$varHIH))[-1]
      }
      or <- exp(co)
      li <- exp(co - z * se)
      ls <- exp(co + z * se)
      r <- cbind(or, li, ls)

      dimnames(r) <- list(names(co), c(lab, paste(level * 100, "%", sep = ""), "C.I."))

      n <- r

      dd <- dim(n)
      n[n > 999.99] <- Inf
      a <- formatC(n, d, len, format = "f")

      dim(a) <- dd
      if (length(dd) == 1) {
        dd <- c(1, dd)
        dim(a) <- dd
        lab <- " "
      } else {
        lab <- dimnames(n)[[1]]
      }

      mx <- max(nchar(lab)) + 1

      cat("Recurrences:\n")
      cat("------------- \n")
      cat(paste(rep(" ", mx), collapse = ""), paste("   ", dimnames(n)[[2]]), "\n")
      for (i in 1:x$nvarnotdep[1])
      {
        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), collapse = "")
        cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }


      cat("\n")
      cat("Terminal event:\n")
      cat("--------------- \n")
      cat(paste(rep(" ", mx), collapse = ""), paste("   ", dimnames(n)[[2]]), "\n")
      for (i in (x$nvarnotdep[1] + 1):dd[1])
      {
        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), collapse = "")
        cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }
    } else {
      x <- object
      if (!inherits(x, "jointPenal")) {
        stop("Object must be of class 'frailtyPenal'")
      }

      z <- abs(qnorm((1 - level) / 2))
      co <- x$coef
      indic_alpha <- x$indic_alpha
      if (indic_alpha == 1 || x$joint.clust == 2) {
        se <- sqrt(diag(x$varHIH))[-c(1, 2)]
      } else {
        se <- sqrt(diag(x$varHIH))[-1]
      }
      or <- co
      li <- co - z * se
      ls <- co + z * se
      r <- cbind(or, li, ls)

      dimnames(r) <- list(names(co), c("coefs", paste(level * 100, "%", sep = ""), "C.I."))

      n <- r

      dd <- dim(n)
      n[n > 999.99] <- Inf
      a <- formatC(n, d, len, format = "f")

      dim(a) <- dd
      if (length(dd) == 1) {
        dd <- c(1, dd)
        dim(a) <- dd
        lab <- " "
      } else {
        lab <- dimnames(n)[[1]]
      }

      mx <- max(nchar(lab)) + 1

      cat("Recurrences:\n")
      cat("------------- \n")
      cat(
        paste(rep(" ", mx - 2), collapse = ""),
        paste("  ", dimnames(n)[[2]][1]),
        paste("", dimnames(n)[[2]][2]),
        paste("", dimnames(n)[[2]][3]),
        "\n"
      )
      # cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
      for (i in 1:x$nvarnotdep[1])
      {
        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), collapse = "")
        cat(lab[i], " ", a[i, 1], "  (", a[i, 2], ",", a[i, 3], " )", sep = "")
        cat("\n")
        # cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }


      cat("\n")
      cat("Terminal event:\n")
      cat("--------------- \n")
      cat(
        paste(rep(" ", mx - 2), collapse = ""),
        paste("  ", dimnames(n)[[2]][1]),
        paste("", dimnames(n)[[2]][2]),
        paste("", dimnames(n)[[2]][3]),
        "\n"
      )
      # cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
      for (i in (x$nvarnotdep[1] + 1):dd[1])
      {
        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), collapse = "")
        cat(lab[i], " ", a[i, 1], "  (", a[i, 2], ",", a[i, 3], " )", sep = "")
        cat("\n")
        # cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }
    }


    cat("\n")

    fmt_pval <- function(p) {
      if (p < 0.001) {
        return("< 0.001")
      }
      return(sprintf("%.4f", p))
    }

    # Check if it is a General Joint Frailty Model (GJFM)
    if (!is.null(x$joint.clust) && x$joint.clust == 2) {
      # theta (u = association rec.-death)
      thetaParam <- x$theta
      thetaSE <- sqrt(((2 * (x$theta^0.5))^2) * diag(x$varHIH)[1]) # delta method
      pvalTheta <- 1 - pnorm(thetaParam / thetaSE)

      if (pvalTheta < 0.05) {
        msg_theta <- "The shared frailty variance theta is significantly different from zero"
        imp_theta <- "a significant association exists between the recurrent events and the terminal event."
      } else {
        msg_theta <- "The shared frailty variance theta is not significantly different from zero"
        imp_theta <- "the recurrent event process and the terminal event process are considered independent."
      }

      thetaInterpret <- sprintf(
        "%s (theta = %.4f, p = %s). Consequently, %s",
        msg_theta, thetaParam, fmt_pval(pvalTheta), imp_theta
      )

      # still GJFM (eta = intra-subject correlation)
      etaParam <- x$eta
      etaSE <- sqrt(((2 * (x$eta^0.5))^2) * diag(x$varHIH)[2]) # delta method
      pvalEta <- 1 - pnorm(etaParam / etaSE)

      if (pvalEta < 0.05) {
        msg_eta <- "The specific frailty variance eta is significantly different from zero"
        imp_eta <- "there is significant unobserved heterogeneity (intra-subject correlation) regarding the recurrent events."
      } else {
        msg_eta <- "The specific frailty variance eta is not significantly different from zero"
        imp_eta <- "there is no significant unobserved heterogeneity regarding the recurrent events."
      }

      etaInterpret <- sprintf(
        "%s (eta = %.4f, p = %s). This indicates that %s",
        msg_eta, etaParam, fmt_pval(pvalEta), imp_eta
      )

      cat(strwrap(thetaInterpret), sep = "\n")
      cat("\n")
      cat(strwrap(etaInterpret), sep = "\n")
      cat("\n")
    } else {
      # Joint Frailty Model (JFM)
      if (x$logNormal == 0) {
        # Gamma Distribution
        frailtyParam <- x$theta
        # Delta method : Var(theta) ~ (2*sqrt(theta))^2 * Var(sqrt(theta))
        frailtySE <- sqrt(((2 * (x$theta^0.5))^2) * diag(x$varHIH)[1])
        label_var <- "theta"
      } else {
        # Log-normal
        frailtyParam <- x$sigma2
        frailtySE <- sqrt(((2 * (x$sigma2^0.5))^2) * diag(x$varHIH)[1])
        label_var <- "sigma2"
      }

      pvalFrailty <- 1 - pnorm(frailtyParam / frailtySE)

      if (pvalFrailty < 0.05) {
        msg_het <- "There is significant unobserved heterogeneity among subjects"
      } else {
        msg_het <- "There is no significant unobserved heterogeneity among subjects"
      }

      thetaInterpret <- sprintf(
        "%s (frailty variance %s = %.4f, p = %s).",
        msg_het, label_var, frailtyParam, fmt_pval(pvalFrailty)
      )

      alphaInterpret <- ""

      if (x$joint.clust != 2) {
        if (x$indic_alpha == 1) {
          # alpha is estimated
          is_signif <- x$alpha_p.value < 0.05
          if (is_signif) {
            msg_alpha <- "The association parameter alpha is significantly different from zero"

            # Only interpret direction if significant
            if (x$alpha > 0) {
              dir_alpha <- "a higher risk of recurrent events is associated with a higher risk of terminal event (alpha > 0)."
            } else {
              dir_alpha <- "a higher risk of recurrent events is associated with a lower risk of terminal event (alpha < 0)."
            }
          } else {
            msg_alpha <- "The association parameter alpha is not significantly different from zero"
            dir_alpha <- "the two processes are considered independent."
          }

          alphaInterpret <- sprintf(
            "%s (alpha = %.4f, p = %s). \nHence, %s",
            msg_alpha, x$alpha, fmt_pval(x$alpha_p.value), dir_alpha
          )
        } else {
          # alpha is fixed (with argument Alpha = "None")
          alphaInterpret <- sprintf(
            "The association parameter was fixed (alpha = 1). A higher risk of recurrent events is associated with a higher risk of terminal event."
          )
        }
      }

      cat(strwrap(thetaInterpret), sep = "\n")
      cat("\n")
      cat(strwrap(alphaInterpret), sep = "\n")
      cat("\n")
    }
  }
