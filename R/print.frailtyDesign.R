#' Print a short table of a 'frailtyDesign' result.
#'
#' @usage \method{print}{frailtyDesign}(x, digits = 2, ...)
#'
#' @param x an object of class 'frailtyDesign' (output from one of the
#' *.power or *.ssize functions).
#' @param digits number of decimals to print for numeric fields. Default is 2.
#' @param \dots other unused arguments.
#'
#' @examples
#' est.ex <- SFM.power(
#'   Groups = 400, ni = 3, ni.type = "max", FUP = 6, Acc.Dur = 0.5, median.H0 = 1.5,
#'   beta.HA = log(0.7), theta = 0.5, cens.par = c(3, 10), cens.type = "Unif", data.type = "rec_event"
#' )
#'
#' print(est.ex)
#'
#' @keywords methods
#' @seealso \code{\link{frailtyDesign}}
#' @export
#' @md
print.frailtyDesign <- function(x, digits = 2, ...) {
  if (!inherits(x, "frailtyDesign")) {
    stop("x must inherit from class 'frailtyDesign'.")
  }

  # --- Helper Functions ---

  # Continued fraction algorithm for obtaining the simplest fraction
  getSimplestFraction <- function(x, tol = 1e-8, maxDenom = 2000) {
    if (x == 0) return(c(numerator = 0, denominator = 1))
    a <- floor(x)
    h1 <- a; h2 <- 1
    k1 <- 1; k2 <- 0
    y <- x - a
    while (y > tol && h1 < maxDenom) {
      y <- 1 / y
      a <- floor(y)
      h <- a * h1 + h2
      k <- a * k1 + k2
      h2 <- h1; k2 <- k1
      h1 <- h; k1 <- k
      y <- y - a
    }
    c(numerator = h1, denominator = k1)
  }

  # Compute the allocation ratio as a simplified fraction
  getAllocRatio <- function(ratio = NULL) {
    if (is.null(ratio) || is.na(ratio)) ratio <- 1
    if (ratio >= 1) {
      frac <- getSimplestFraction(ratio)
      list(trt = frac["numerator"], ctrl = frac["denominator"])
    } else {
      frac <- getSimplestFraction(1 / ratio)
      list(trt = frac["denominator"], ctrl = frac["numerator"])
    }
  }

  # Numeric formatting
  fmtNum <- function(val, d = digits) formatC(val, digits = d, format = "f")

  # Build a label/value pair for the 'ni' parameter
  buildNiLine <- function(model, dType, niVal, niType) {
    if (anyNA(niVal)) return(NULL)

    interpretRecEvents <- function(val, valType) {
      if (is.null(valType) || is.na(valType)) {
        c(label = "Rec. events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("m", "max", "maximum")) {
        c(label = "Maximum number of rec. events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("pois", "poisson", "p")) {
        c(label = "Mean number of rec. events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("u", "unif", "uniform") && length(val) == 2) {
        c(label = "Number of rec. events per subject range",
          value = paste(val[1], "to", val[2]))
      } else {
        c(label = "Rec. events per subject (unknown type)", value = as.character(val))
      }
    }

    if (model == "SFM") {
      if (dType == "grouped") return(c(label = "Subjects per group", value = as.character(niVal)))
      if (dType == "rec_event") return(interpretRecEvents(niVal, niType))
    }
    if (model == "NFM") {
      if (dType == "grouped") return(c(label = "Subgroups per group", value = as.character(niVal)))
      if (dType == "rec_event2") return(c(label = "Number of event types per subject", value = as.character(niVal)))
      if (dType == "rec_event1") return(c(label = "Subjects per group", value = as.character(niVal)))
    }
    if (model %in% c("JFM", "GJFM")) return(interpretRecEvents(niVal, niType))

    c(label = "ni", value = as.character(niVal))
  }

  # Build a label/value pair for the 'kij' parameter
  buildKijLine <- function(model, dType, kijVal, kijType) {
    if (anyNA(kijVal)) return(NULL)

    interpretRecEvents <- function(val, valType, isRecEvent2 = FALSE) {
      baseLabel <- if (isRecEvent2) "rec. events for each type of event" else "rec. events per subject"
      if (is.null(valType) || is.na(valType)) {
        c(label = paste("Number of", baseLabel), value = as.character(val))
      } else if (tolower(valType) %in% c("m", "max", "maximum")) {
        c(label = paste("Maximum number of", baseLabel), value = as.character(val))
      } else if (tolower(valType) %in% c("pois", "poisson", "p")) {
        c(label = paste("Mean number of", baseLabel), value = as.character(val))
      } else if (tolower(valType) %in% c("u", "unif", "uniform") && length(val) == 2) {
        c(label = paste("Number of", baseLabel, "ranging from"),
          value = paste(val[1], "to", val[2]))
      } else {
        c(label = paste("Number of", baseLabel, "(unknown type)"), value = as.character(val))
      }
    }

    if (model == "SFM") {
      if (dType == "grouped") return(c(label = "Number of subjects per subgroup", value = as.character(kijVal)))
      if (dType == "rec_event") return(interpretRecEvents(kijVal, kijType))
    }
    if (model == "NFM") {
      if (dType == "grouped") return(c(label = "Number of subjects per subgroup", value = as.character(kijVal)))
      if (dType == "rec_event1") return(interpretRecEvents(kijVal, kijType, isRecEvent2 = FALSE))
      if (dType == "rec_event2") return(interpretRecEvents(kijVal, kijType, isRecEvent2 = TRUE))
    }
    if (model %in% c("JFM", "GJFM")) return(interpretRecEvents(kijVal, kijType))

    c(label = "kij", value = as.character(kijVal))
  }

  # Helper to print a formatted line (label and value)
  prLine <- function(label, value, widthLabel = 52) {
    cat(sprintf("%-*s: %s\n", widthLabel, label, value))
  }

  # --- Extract Key Information ---
  method    <- x$method
  modelName <- if (!is.null(x$model)) x$model else "Unknown model"
  dType     <- if (!is.null(x$data.type)) x$data.type else NA
  ratioVal  <- if (!is.null(x$ratio)) x$ratio else 1

  # Groups, ni, kij, total subjects
  groupCount <- if (!is.null(x$Groups)) x$Groups else NA
  niVal      <- if (!is.null(x$ni)) x$ni else NA
  niType     <- if (!is.null(x$ni.type)) x$ni.type else NA
  kijVal     <- if (!is.null(x$kij)) x$kij else NA
  kijType    <- if (!is.null(x$kij.type)) x$kij.type else NA
  totalS     <- if (!is.null(x$Npts)) x$Npts else if (x$data.type != "grouped" && x$model == "SFM" && !is.null(x$Groups)) x$Groups else NA

  # Decide scenario based on model and data type
  decideScenario <- function(m, dt) {
    if (m == "SFM") {
      if (dt == "grouped") return("groupBased")
      if (dt == "rec_event") return("subjectBased")
    } else if (m == "NFM") {
      if (dt %in% c("grouped", "rec_event1")) return("groupBased")
      if (dt == "rec_event2") return("nfmRecEvent2")
    } else if (m %in% c("JFM", "GJFM")) {
      return("subjectBased")
    }
    "unknown"
  }
  scenario <- decideScenario(modelName, dType)

  # --- Compute Sample Size Details ---
  mainSampleSize <- NA
  nSubjTrt       <- NA
  nSubjCtrl      <- NA
  totalSubj      <- NA

  if (scenario == "groupBased") {
    if (!is.na(groupCount)) {
      mainSampleSize <- groupCount
      # For grouped designs the total subjects is computed from groups and (ni, and possibly kij)
      if (modelName == "NFM" && dType == "grouped" && !is.null(x$kij)) {
        totalSubj <- groupCount * niVal * kijVal
      } else {
        totalSubj <- groupCount * niVal
      }
      nSubjTrt <- round(totalSubj * ratioVal / (1 + ratioVal))
      nSubjCtrl <- totalSubj - nSubjTrt
    }
  } else if (scenario == "subjectBased") {
    if (!is.na(totalS)) {
      mainSampleSize <- totalS
      totalSubj <- totalS
      nSubjTrt <- round(totalS * ratioVal / (1 + ratioVal))
      nSubjCtrl <- totalS - nSubjTrt
    }
  } else if (scenario == "nfmRecEvent2") {
    if (!is.na(groupCount)) {
      mainSampleSize <- groupCount
      totalSubj <- groupCount
      nSubjTrt <- round(groupCount * ratioVal / (1 + ratioVal))
      nSubjCtrl <- groupCount - nSubjTrt
    }
  }

  # --- Build Table Lines for Recurrent Event Parameters ---
  niLine  <- buildNiLine(modelName, dType, niVal, niType)
  kijLine <- buildKijLine(modelName, dType, kijVal, kijType)

  # --- Title and Short Description ---
  cat("\n")
  if (method == "power") {
    cat(sprintf("Power Estimation Using a %s\n", modelName))
  } else if (method == "ssize") {
    cat(sprintf("Sample Size Estimation Using a %s\n", modelName))
  }

  modelDescription <- function(model, dType) {
    if (model == "NFM" && !is.null(dType)) {
      if (dType == "grouped") return("Two levels of clustering")
      if (dType == "rec_event1") return("Recurrent event among clusters")
      if (dType == "rec_event2") return("Multi-type recurrent events")
    }
    if (model == "SFM" && !is.null(dType)) {
      if (dType == "grouped") return("Groups and subjects per group")
      if (dType == "rec_event") return("Recurrent events")
    }
    NULL
  }
  shortDesc <- modelDescription(modelName, dType)
  if (!is.null(shortDesc)) cat(shortDesc, "\n")
  cat(strrep("-", 60), "\n")

  ## A) TOP-LEVEL: Print estimated power, group-level details and subject counts

  if (method == "power" && !is.null(x$estimated.power)) {
    prLine("Estimated power", paste0(fmtNum(100 * x$estimated.power), "%"))
  }

  # Print total number of groups if applicable
  if (!is.na(groupCount) && !dType %in% c("rec_event2", "rec_event")) {
    prLine("Total number of groups", as.character(groupCount))
  }

  # Print recurrent event details (ni then kij) if applicable
  if (!is.null(niLine)) prLine(niLine["label"], niLine["value"])
  if (!is.null(kijLine)) prLine(kijLine["label"], kijLine["value"])

  # Then print total number of subjects and allocation breakdown
  if (!is.na(totalSubj)) {
    prLine("Total number of subjects", as.character(totalSubj))
    if (!is.na(nSubjTrt)) prLine("  Experimental", as.character(nSubjTrt))
    if (!is.na(nSubjCtrl)) prLine("  Control", as.character(nSubjCtrl))
  }

  cat("\n")

  ## B) NUMBER OF EVENTS
  if (!is.null(x$number.events) && length(x$number.events) == 2) {
    cat("Expected total number of events\n")
    prLine("  Under the null (H0)", sprintf("%d", ceiling(x$number.events[1])))
    prLine("  Under the alternative (HA)", sprintf("%d", ceiling(x$number.events[2])))
    cat("\n")
  }
  if (!is.null(x$number.events.rec) && length(x$number.events.rec) == 2) {
    cat("Expected total number of rec. events\n")
    prLine("  Under the null (H0)", sprintf("%d", ceiling(x$number.events.rec[1])))
    prLine("  Under the alternative (HA)", sprintf("%d", ceiling(x$number.events.rec[2])))
    cat("\n")
  }
  if (!is.null(x$number.events.D) && length(x$number.events.D) == 2) {
    cat("Expected total number of terminal events\n")
    prLine("  Under the null (H0)", sprintf("%d", ceiling(x$number.events.D[1])))
    prLine("  Under the alternative (HA)", sprintf("%d", ceiling(x$number.events.D[2])))
    cat("\n")
  }

  ## C) PARAMETERS
  cat("With:\n")

  if (!is.null(x$alpha)) {
    prLine("  Significance level", paste0(fmtNum(100 * x$alpha), "%"))
  }
  if (!is.null(x$testType)) {
    prLine("  Test type", x$testType)
  }
  if (!is.null(x$Acc.Dur)) {
    prLine("  Accrual time", x$Acc.Dur)
  }
  if (!is.null(x$FUP)) {
    prLine("  Follow-up time", x$FUP)
  }
  if (!is.null(x$FUP.type)) {
    prLine("  Subjects are followed-up", ifelse(x$FUP.type == "fixed", "for a fixed period", "until the end of study"))
  }
  if (!is.null(x$tested.structure)) {
    testedStr <- ifelse(tolower(x$tested.structure) == "joint",
                        "both recurrent and terminal events",
                        ifelse(tolower(x$tested.structure) == "betartest",
                               "the recurrent events only",
                               "the terminal events only"))
    prLine("  Tested structure", testedStr)
  }
  if (!is.null(x$ratio)) {
    allocRatio <- getAllocRatio(x$ratio)
    prLine("  Allocation ratio (exp:ctrl)", paste0(as.numeric(allocRatio$trt), ":", as.numeric(allocRatio$ctrl)))
  }
  if (!is.null(x$timescale) && (dType != "grouped" | is.na(dType))) {
    prLine("  Timescale", x$timescale)
  }
  
  # Added frailty variance parameters:
  if (!is.null(x$theta)) {
    prLine("  Main gamma-frailty variance", fmtNum(x$theta))
  }
  if (modelName %in% c("NFM", "GJFM") && !is.null(x$eta)) {
    if (modelName == "NFM") {
      prLine("  Second-level nesting variance", fmtNum(x$eta))
    } else if (modelName == "GJFM") {
      prLine("  Variance of inter-reccurence dependance", fmtNum(x$eta))
    }
  }
  
  if (!is.null(x$samplesMC)) {
    prLine("  Monte Carlo samples", x$samplesMC)
  }

  cat("\n")

  ## D) NOTE: Total number of subjects
  if (method %in% c("ssize", "power")) {
    if (modelName == "SFM" && dType == "grouped" && !is.null(x$Groups) && !is.null(x$ni)) {
      cat(sprintf("Note: This corresponds to a total of %d subjects (%d groups x %d subjects per group).\n",
                  x$Groups * x$ni, x$Groups, x$ni))
    } else if (modelName == "SFM" && dType == "rec_event" && !is.null(x$Npts)) {
      cat(sprintf("Note: This corresponds to a total of %d subjects (as provided by Npts for SFM rec_event).\n", x$Npts))
    } else if (modelName == "NFM" && dType == "grouped" && !is.null(x$Groups) && !is.null(x$ni) && !is.null(x$kij)) {
      cat(sprintf("Note: This leads to a total of %d subjects (%d groups x %d subgroups per group x %d subjects per subgroup).\n",
                  x$Groups * x$ni * x$kij, x$Groups, x$ni, x$kij))
    } else if (modelName == "NFM" && dType == "rec_event1" && !is.null(x$Groups) && !is.null(x$ni)) {
      cat(sprintf("Note: This leads to a total of %d subjects (%d groups x %d subjects per group).\n",
                  x$Groups * x$ni, x$Groups, x$ni))
    }
  }

  cat("\n")
  invisible(x)
}
