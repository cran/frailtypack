#' Summarize a 'frailtyDesign' object.
#'
#' @usage \method{summary}{frailtyDesign}(object, digits = 2, ...)
#'
#' @param object an object of class 'frailtyDesign' (output from one of the
#'   *.power or *.ssize functions).
#' @param digits number of decimals to print for numeric fields. Default is 2.
#' @param \dots other unused arguments.
#'
#' @examples
#' est.ex <- SFM.power(
#'   Groups = 400, ni = 3, ni.type = "max", FUP = 6, Acc.Dur = 0.5, median.H0 = 1.5,
#'   beta.HA = log(0.7), theta = 0.5, cens.par = c(3, 10), cens.type = "Unif", data.type = "rec_event"
#' )
#'
#' summary(est.ex)
#'
#' @keywords methods
#' @seealso \code{\link{frailtyDesign}}
#' @export
#' @md
summary.frailtyDesign <- function(object, digits = 2, ...) {
  x <- object
  if (!inherits(x, "frailtyDesign")) {
    stop("x must inherit from class 'frailtyDesign'.")
  }

  # --- Helper Functions ---

  # continued fraction algorithm
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

  # Compute the allocation ratio (treatment:control) as a simplified fraction
  getAllocRatio <- function(ratio = NULL) {
    if (is.null(ratio) || is.na(ratio)) ratio <- 1
    if (ratio >= 1) {
      frac <- getSimplestFraction(ratio)
      treatment <- frac["numerator"]
      control <- frac["denominator"]
    } else {
      invRatio <- 1 / ratio
      frac <- getSimplestFraction(invRatio)
      treatment <- frac["denominator"]
      control <- frac["numerator"]
    }
    list(trt = treatment, ctrl = control)
  }

  formatHR <- function(hr) {
    if (is.null(hr) || is.na(hr)) return(NULL)
    if (abs(hr - 1) < 1e-12) {
      "1.0 (lack of association)"
    } else if (hr < 1) {
      sprintf("%.2f (%.1f%% risk reduction)", hr, 100 * (1 - hr))
    } else {
      sprintf("%.2f (%.1f%% risk increase)", hr, 100 * (hr - 1))
    }
  }

  abbreviateRec <- function(label) {
    # Replace "recurrent" with "rec." (case-insensitive) in the given label
    gsub("recurrent", "rec.", label, ignore.case = TRUE)
  }

  # Build a label/value pair for the 'ni' parameter
  buildNiLine <- function(model, dType, niVal, niType) {
    if (anyNA(niVal)) return(NULL)

    interpretRecEvents <- function(val, valType) {
      if (is.null(valType) || is.na(valType)) {
        c(label = "Recurrent events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("m", "max", "maximum")) {
        c(label = "Maximum number of rec. events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("pois", "poisson", "p")) {
        c(label = "Mean number of rec. events per subject", value = as.character(val))
      } else if (tolower(valType) %in% c("u", "unif", "uniform") && length(val) == 2) {
        c(label = "Number of rec. events per subject ranging from",
          value = paste(val[1], "to", val[2]))
      } else {
        c(label = "Recurrent events per subject (unknown type)", value = as.character(val))
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

    c(label = "niVal", value = as.character(niVal))
  }

  # Build a label/value pair for the 'kij' parameter
  buildKijLine <- function(model, dType, kijVal, kijType) {
    if (anyNA(kijVal)) return(NULL)

    interpretRecEvents <- function(val, valType, isRecEvent2 = FALSE) {
      baseLabel <- if (isRecEvent2) "recurrent events for each type of event" else "recurrent events per subject"
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

    c(label = "kijVal", value = as.character(kijVal))
  }

  # Print a formatted line (label and value)
  prLine <- function(label, value, widthLabel = 52) {
    cat(sprintf("%-*s: %s\n", widthLabel, label, value))
  }

  # Format percentages neatly
  .formatPerc <- function(x, maxDecimals = 2) {
    xRounded <- round(x, maxDecimals)
    if (abs(xRounded - floor(xRounded)) < 1e-12) {
      paste0(as.integer(xRounded), "%")
    } else {
      paste0(xRounded, "%")
    }
  }

  # --- Extract Key Information ---
  method   <- x$method        # "power" or "ssize"
  model    <- x$model         # "SFM", "NFM", "JFM", or "GJFM"
  alpha    <- if (!is.null(x$alpha)) x$alpha else NA
  tsc      <- if (!is.null(x$timescale)) x$timescale else NA
  ratioVal <- if (!is.null(x$ratio)) x$ratio else 1
  dType    <- if (!is.null(x$data.type)) x$data.type else NA
  testType <- if (!is.null(x$testType)) x$testType else "2-sided"
  sideText <- if (tolower(testType) == "1-sided") "1-sided" else "2-sided"

  testedStr <- if (!is.null(x$tested.structure)) {
    ifelse(tolower(x$tested.structure) == "joint",
           "both recurrent and terminal events",
           ifelse(tolower(x$tested.structure) == "betartest",
                  "the recurrent events only",
                  "the terminal events only"))
  } else {
    NA
  }

  # Hazard ratios
  HRrecur.H0 <- if (!is.null(x$HR.H0)) x$HR.H0 else x$HR.R0
  HRrecur.HA <- if (!is.null(x$HR.HA)) x$HR.HA else x$HR.RA
  HRterm.H0  <- if (!is.null(x$HR.D0)) x$HR.D0 else NA
  HRterm.HA  <- if (!is.null(x$HR.DA)) x$HR.DA else NA

  # Events count (H0 & HA)
  nEv.recH0 <- if (!is.null(x$number.events)) x$number.events[1] else NA
  nEv.recHA <- if (!is.null(x$number.events)) x$number.events[2] else NA
  if (!is.null(x$number.events.rec)) {
    nEv.recH0 <- x$number.events.rec[1]
    nEv.recHA <- x$number.events.rec[2]
  }
  nEv.dthH0 <- if (!is.null(x$number.events.D)) x$number.events.D[1] else NA
  nEv.dthHA <- if (!is.null(x$number.events.D)) x$number.events.D[2] else NA

  # Groups, ni, kij, total subjects
  groupCount <- if (!is.null(x$Groups)) x$Groups else NA
  niVal      <- if (!is.null(x$ni)) x$ni else NA
  niType     <- if (!is.null(x$ni.type)) x$ni.type else NA
  kijVal     <- if (!is.null(x$kij)) x$kij else NA
  kijType    <- if (!is.null(x$kij.type)) x$kij.type else NA

  totalS <- if (!is.null(x$Npts)) x$Npts else NA
  if (dType != "grouped" && model == "SFM" && !is.na(groupCount) && !anyNA(niVal))
    totalS <- groupCount

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
  scenario <- decideScenario(model, dType)

  # --- Compute Sample Size Details ---
  mainSampleSize <- NA
  nSubjTrt       <- NA
  nSubjCtrl      <- NA
  nGroupsTrt     <- NA
  nGroupsCtrl    <- NA
  totalSubj      <- NA

  if (scenario == "groupBased") {
    if (!is.na(groupCount)) {
      mainSampleSize <- groupCount
      if (!anyNA(niVal) && !anyNA(kijVal) && dType == "grouped") {
        totalSubj <- groupCount * niVal * kijVal
      } else if (!anyNA(niVal)) {
        totalSubj <- groupCount * niVal
      }
      nGroupsTrt <- round(groupCount * ratioVal / (1 + ratioVal))
      nGroupsCtrl <- groupCount - nGroupsTrt
      if (!is.na(totalSubj)) {
        ratioSplit <- nGroupsTrt / groupCount
        nSubjTrt <- round(totalSubj * ratioSplit)
        nSubjCtrl <- totalSubj - nSubjTrt
      }
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
      nGroupsTrt <- round(groupCount * ratioVal / (1 + ratioVal))
      nGroupsCtrl <- groupCount - nGroupsTrt
      nSubjTrt <- nGroupsTrt
      nSubjCtrl <- nGroupsCtrl
    }
  }

  # --- Build Table Lines ---
  niLine  <- buildNiLine(model, dType, niVal, niType)
  kijLine <- buildKijLine(model, dType, kijVal, kijType)

  # --- Title and Short Description ---
  if (method == "power") {
    cat(sprintf("Power Estimation Using a %s\n", model))
  } else {
    cat(sprintf("Sample Size Estimation Using a %s\n", model))
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
  shortDesc <- modelDescription(model, dType)
  if (!is.null(shortDesc)) cat(shortDesc, "\n")
  cat(strrep("-", 60), "\n")

  ## --- Table Output ---
  if (method == "power" && !is.null(x$estimated.power)) {
    prLine("Estimated power", .formatPerc(100 * x$estimated.power))
  }

  if (!is.na(groupCount) && !(dType %in% c("rec_event2", "rec_event"))) {
    prLine("Total number of groups", as.character(groupCount))
  }

  if (!is.null(niLine)) prLine(abbreviateRec(niLine["label"]), niLine["value"])
  if (!is.null(kijLine)) prLine(abbreviateRec(kijLine["label"]), kijLine["value"])

  if (!is.na(totalSubj)) {
    prLine("Total number of subjects", as.character(totalSubj))
    if (!is.na(nSubjTrt)) prLine("  Experimental", as.character(nSubjTrt))
    if (!is.na(nSubjCtrl)) prLine("  Control", as.character(nSubjCtrl))
  }

  cat("\n")

  # --- Build Summary Sentences ---

  makeRecEventParenthetical <- function(lbl, val) {
    lblLow <- tolower(lbl)
    if (grepl("maximum number of (?:rec\\.?|recurrent) events per subject", lblLow, perl = TRUE)) {
      sprintf("and a maximum number of %s recurrent events per subject", val)
    } else if (grepl("mean number of (?:rec\\.?|recurrent) events per subject", lblLow, perl = TRUE)) {
      sprintf("and a mean number of %s recurrent events per subject", val)
    } else if (grepl("number of (?:rec\\.?|recurrent) events per subject ranging from", lblLow, perl = TRUE)) {
      sprintf("and a number of recurrent events per subject ranging from %s", val)
    } else if (grepl("(?:rec\\.?|recurrent) events per subject", lblLow, perl = TRUE)) {
      sprintf("and %s %s", lbl, val)
    } else if (grepl("(?:rec\\.?|recurrent) events for each type of event", lblLow, perl = TRUE)) {
      if (grepl("maximum number", lblLow)) {
        sprintf("and a maximum number of %s recurrent events for each type of event", val)
      } else if (grepl("mean number", lblLow)) {
        sprintf("and a mean number of %s recurrent events for each type of event", val)
      } else if (grepl("ranging from", lblLow)) {
        sprintf("and a number of recurrent events for each type of event ranging from %s", val)
      } else {
        sprintf("and %s %s", lbl, val)
      }
    } else {
      ""
    }
  }

  # Modify the check below to capture both "rec." and "recurrent"
  recEventClause <- ""
  if (dType %in% c("rec_event", "rec_event1", "rec_event2") || model %in% c("JFM", "GJFM")) {
    if (!is.null(niLine) && grepl("(rec\\.?|recurrent) event", tolower(niLine["label"]), perl = TRUE)) {
      recEventClause <- makeRecEventParenthetical(niLine["label"], niLine["value"])
    }
    if (!is.null(kijLine) && grepl("(rec\\.?|recurrent) event", tolower(kijLine["label"]), perl = TRUE)) {
      part <- makeRecEventParenthetical(kijLine["label"], kijLine["value"])
      if (nzchar(part)) {
        recEventClause <- if (nzchar(recEventClause))
          paste0(recEventClause, " ", sub("^\\(", "(", part))
        else part
      }
    }
  }

  designLabel <- ""
  if (model == "SFM") {
    if (dType == "grouped" && !is.na(groupCount) && !anyNA(niVal)) {
      designLabel <- sprintf("%d groups, each with %d subjects, for a total of %d subjects",
                             groupCount, niVal, groupCount * niVal)
    } else if (dType == "rec_event" && !is.na(totalS)) {
      designLabel <- sprintf("%d subjects", totalS)
    }
  } else if (model == "NFM") {
    if (dType == "grouped" && !is.na(groupCount) && !anyNA(niVal) && !anyNA(kijVal)) {
      designLabel <- sprintf("%d groups, each with %d subgroups of size %d, for a total of %d subjects",
                             groupCount, niVal, kijVal, groupCount * niVal * kijVal)
    } else if (dType == "rec_event1" && !is.na(groupCount) && !anyNA(niVal)) {
      designLabel <- sprintf("%d groups, each with %d subjects, for a total of %d subjects",
                             groupCount, niVal, groupCount * niVal)
    } else if (dType == "rec_event2" && !is.na(groupCount) && !anyNA(niVal)) {
      designLabel <- sprintf("%d subjects, each with %d distinct type of events",
                             groupCount, niVal)
    }
  } else if (model %in% c("JFM", "GJFM")) {
    if (!is.na(totalS)) designLabel <- sprintf("%d subjects", totalS)
  }
  if (designLabel == "") designLabel <- "the specified design"

  if (nzchar(recEventClause)) {
    if (model == "NFM" && dType == "rec_event2") {
      designLabel <- paste0(designLabel, ", ", recEventClause)
    } else if (grepl("(\\d+) subjects", designLabel)) {
      designLabel <- sub("([0-9]+ subjects)", paste0("\\1 ", recEventClause), designLabel)
    } else if (grepl("(\\d+) groups", designLabel)) {
      designLabel <- sub("([0-9]+ groups)", paste0("\\1 ", recEventClause), designLabel)
    } else {
      designLabel <- paste(designLabel, recEventClause)
    }
  }

  modelClause <- sprintf("Using a %s model", model)
  if (!is.null(testedStr) && !is.na(testedStr) && nzchar(testedStr)) {
    modelClause <- sprintf("%s, testing the treatment effect for %s", modelClause, testedStr)
  }

  opener <- ""
  if (method == "power") {
    if (!is.null(x$estimated.power)) {
      opener <- sprintf("%s, a study with %s, will yield an estimated power of %s.",
                        modelClause, designLabel, .formatPerc(100 * x$estimated.power))
    } else {
      opener <- sprintf("%s, power estimation is based on %s (no numerical power available).",
                        modelClause, designLabel)
    }
  } else if (method == "ssize") {
    tgtP <- if (!is.null(x$target.power)) x$target.power else NA
    if (!is.na(tgtP)) {
      opener <- sprintf("%s, a study with %s, is required to achieve a target power of %s.",
                        modelClause, designLabel, .formatPerc(100 * tgtP))
    } else {
      opener <- sprintf("%s, the required sample size was determined based on %s.",
                        modelClause, designLabel)
    }
  } else {
    opener <- sprintf("%s with %s.", modelClause, designLabel)
  }

  # Assemble additional design assumptions
  bits <- c()
  if (!is.null(tsc) && !is.na(tsc)) {
    if ((model %in% c("SFM", "NFM") && !is.null(dType) && !is.na(dType) && dType != "grouped") ||
         model %in% c("JFM", "GJFM"))
      bits <- c(bits, sprintf("a '%s' timescale", tsc))
  }
  if (!is.na(alpha))
    bits <- c(bits, sprintf("a type-I error of %s (%s)", .formatPerc(100 * alpha), sideText))

  accrualTime <- if (!is.null(x$Acc.Dur)) x$Acc.Dur else NA
  followUpTime <- if (!is.null(x$FUP)) x$FUP else NA
  followUpType <- if (!is.null(x$FUP.type)) ifelse(x$FUP.type == "fixed", "for a fixed period", "until the end of study") else NA

  if (!is.na(followUpTime) && !is.na(followUpType)) {
    bits <- c(bits, sprintf("an accrual time of %s time unit(s) and a follow-up time of %s time unit(s) (follow-up %s)",
                              as.character(accrualTime), as.character(followUpTime), as.character(followUpType)))
  }

  allocRatio <- getAllocRatio(ratioVal)
  if (ratioVal != 1) {
    bits <- c(bits, sprintf("an allocation ratio of %d:%d (experimental:control)", allocRatio$trt, allocRatio$ctrl))
  } else {
    bits <- c(bits, "an equal allocation ratio (1:1)")
  }

  secondPart <- if (length(bits) > 0)
    paste0("This design assumes ", paste(bits, collapse = ", "), ".")
  else ""

  # Hazard ratios and expected events
  hrPart <- ""
  if (model %in% c("SFM", "NFM")) {
    if (!is.null(HRrecur.H0) && !is.null(HRrecur.HA)) {
      hrPart <- sprintf("Under the null hypothesis (H0), the hazard ratio is %s, whereas under the alternative (HA), it is %s.",
                        formatHR(HRrecur.H0), formatHR(HRrecur.HA))
    }
  } else {
    hrPart <- sprintf("Under the null hypothesis (H0), the HR for recurrent event is %s and the HR for terminal event is %s; under the alternative (HA), the HRs are %s and %s, respectively.",
                      formatHR(HRrecur.H0), formatHR(HRterm.H0),
                      formatHR(HRrecur.HA), formatHR(HRterm.HA))
  }

  eventsPart <- ""
  if (model %in% c("SFM", "NFM") && !is.na(nEv.recH0) && !is.na(nEv.recHA)) {
    eventsPart <- sprintf("Based on these assumptions, we anticipate approximately %d events under H0 and %d under HA.",
                          nEv.recH0, nEv.recHA)
  } else if (!is.na(nEv.recH0) && !is.na(nEv.recHA) && !is.na(nEv.dthH0) && !is.na(nEv.dthHA)) {
    eventsPart <- sprintf("Based on these assumptions, we anticipate about %d recurrent events and %d terminal events under H0, versus %d and %d under HA, respectively.",
                          nEv.recH0, nEv.dthH0, nEv.recHA, nEv.dthHA)
  }

  finalReport <- paste(opener, secondPart, hrPart, eventsPart)
  finalReport <- gsub("\\s+", " ", finalReport)  # clean up extra spaces
  cat(finalReport, "\n\n")

  invisible(x)
}

