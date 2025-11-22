#' Breast Cosmesis Data
#' 
#' The often used data set for interval-censored data, described and given in
#' full in Finkelstein and Wolfe (1985). It involves 94 breast cancer patients
#' who were randomized to either radiation therapy with chemotherapy or
#' radiation therapy alone. The outcome is time until the onset of breast
#' retraction which is interval-censored between the last clinic visit before
#' the event was observed and the first visit when the event was observed.
#' Patients without breast retraction were right-censored.
#' 
#' 
#' @name bcos
#' @docType data
#' @usage data(bcos)
#' @format A data frame with 94 observations and 3 variables: \describe{
#' \item{left}{left end point of the breast retraction interval}
#' \item{right}{right end point of the breast retraction interval}
#' \item{treatment}{type of treatment received} }
#' @source Finkelstein, D.M. and Wolfe, R.A. (1985). A semiparametric model for
#' regression analysis of interval-censored failure time data.
#' \emph{Biometrics} \bold{41}, 731-740.
#' @keywords datasets
NULL





#' Follow-up of metastatic colorectal cancer patients: times of new lesions
#' appearance and death
#' 
#' Randomly chosen 150 patients from the follow-up of the FFCD 2000-05
#' multicenter phase III clinical trial originally including 410 patients with
#' metastatic colorectal cancer randomized into two therapeutic strategies:
#' combination and sequential. The dataset contains times of observed
#' appearances of new lesions censored by a terminal event (death or
#' right-censoring) with baseline characteristics (treatment arm, age, WHO
#' performance status and previous resection).
#' 
#' 
#' @name colorectal
#' @docType data
#' @usage data(colorectal)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{time0}{start of interval (0 or previous recurrence time)}
#' \item{time1}{recurrence or censoring time} \item{new.lesions}{Appearance of
#' new lesions status. 0: censsored or no event, 1: new lesions}
#' \item{treatment}{To which treatment arm a patient was allocated? 1:
#' sequential (S); 2: combination (C)} \item{age}{Age at baseline: 1: <50
#' years, 2: 50-69 years, 3: >69 years} \item{who.PS}{WHO performance status at
#' baseline: 1: status 0, 2: status 1, 3: status 2}
#' \item{prev.resection}{Previous resection of the primate tumor?  0: No, 1:
#' Yes} \item{state}{death indicator. 0: alive, 1: dead}
#' \item{gap.time}{interocurrence time or censoring time} }
#' @note We thank the Federation Francophone de Cancerologie Digestive and
#' Gustave Roussy for sharing the data of the FFCD 2000-05 trial supported by
#' an unrestricted Grant from Sanofi.
#' @references M. Ducreux, D. Malka, J. Mendiboure, P.-L. Etienne, P. Texereau,
#' D. Auby, P. Rougier, M. Gasmi, M. Castaing, M. Abbas, P. Michel, D. Gargot,
#' A. Azzedine, C. Lombard- Bohas, P. Geoffroy, B. Denis, J.-P., Pignon,
#' L.,Bedenne, and O.  Bouche (2011). Sequential versus combination
#' chemotherapy for the treatment of advanced colorectal cancer (FFCD 2000-05):
#' an open-label, randomised, phase 3 trial.  \emph{The Lancet Oncology}
#' \bold{12}, 1032-44.
#' @keywords datasets
NULL





#' Follow-up of metastatic colorectal cancer patients : longitudinal
#' measurements of tumor size
#' 
#' Randomly chosen 150 patients from the follow-up of the FFCD 2000-05
#' multicenter phase III clinical trial originally including 410 patients with
#' metastatic colorectal cancer randomized into two therapeutic strategies:
#' combination and sequential. The dataset contains measurements of tumor size
#' (left-censored sums of the longest diameters of target lesions; transformed
#' using Box-Cox) with baseline characteristics(treatment arm, age, WHO
#' performance status and previous resection).
#' 
#' 
#' @name colorectalLongi
#' @docType data
#' @usage data(colorectalLongi)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{year}{time of visit counted in years from baseline}
#' \item{tumor.size}{Individual longitudinal measurement of transformed
#' (Box-Cox with parameter 0.3) sums of the longest diameters, left-censored
#' due to a detection limit (threshold \eqn{s=-3.33}). } \item{treatment}{To
#' which treatment arm a patient was allocated? 1: sequential (S); 2:
#' combination (C)} \item{age}{Age at baseline: 1: <50 years, 2: 50-69 years,
#' 3: >69 years} \item{who.PS}{WHO performance status at baseline: 1: status 0,
#' 2: status 1, 3: status 2} \item{prev.resection}{Previous resection of the
#' primate tumor?  0: No, 1: Yes} }
#' @note We thank the Federation Francophone de Cancerologie Digestive and
#' Gustave Roussy for sharing the data of the FFCD 2000-05 trial supported by
#' an unrestricted Grant from Sanofi.
#' @references Ducreux, M., Malka, D., Mendiboure, J., Etienne, P.-L.,
#' Texereau, P., Auby, D., Rougier, P., Gasmi, M., Castaing, M., Abbas, M.,
#' Michel, P., Gargot, D., Azzedine, A., Lombard- Bohas, C., Geoffroy, P.,
#' Denis, B., Pignon, J.-P., Bedenne, L., and Bouche, O. (2011). Sequential
#' versus combination chemotherapy for the treatment of advanced colorectal
#' cancer (FFCD 2000-05): an open-label, randomised, phase 3 trial.  \emph{The
#' Lancet Oncology} \bold{12}, 1032-44.
#' @keywords datasets
NULL





#' Simulated data as a gathering of clinical trials databases
#' 
#' This contains simulated samples of 100 clusters with 100 subjects in each
#' cluster, like a gathering of clinical trials databases. Two correlated
#' centred gaussian random effects are generated with the same variance fixed
#' at 0.3 and the covariance at -0.2. The regression coefficient \eqn{\beta} is
#' fixed at -0.11. The percentage of right-censored data is around 30 percent
#' which are generated from a uniform distribution on [1,150]. The baseline
#' hazard function is considered as a simple Weibull.
#' 
#' 
#' @name dataAdditive
#' @docType data
#' @usage data(dataAdditive)
#' @format This data frame contains the following columns: \describe{
#' \item{group}{identification variable} \item{t1}{start of interval (=0,
#' because left-truncated data are not allowed)} \item{t2}{end of interval
#' (death or censoring time)} \item{event}{censoring status (0:alive, 1:death,
#' as acensoring indicator} \item{var1}{dichotomous covariate (=0 or 1,as a
#' treatment variable)} \item{var2}{dichotomous covariate (=0 or 1,as a
#' treatment variable)} }
#' @source
#' 
#' V. Rondeau, S. Michiels, B.Liquet, and J.P. Pignon (2008). Investigating
#' trial and treatment heterogeneity in an individual patient data
#' meta-analysis of survival data by mean of the maximum penalized likelihood
#' approach. \emph{Statistics in Medecine}, \bold{27}, 1894-1910.
#' @keywords datasets
NULL





#' Simulated data for two types of recurrent events and a terminal event
#' 
#' This contains a simulated sample of of 800 subjects and 1652 observations.
#' This dataset can be used to illustrate how to fit a joint multivariate
#' frailty model. Two gaussian correlated random effects were generated with
#' mean 0, variances 0.5 and a correlation coefficient equals to 0.5. The
#' coefficients \eqn{\alpha_1} and \eqn{\alpha_2} were fixed to 1. The three
#' baseline hazard functions followed a Weibull distribution and right
#' censoring was fixed at 5.
#' 
#' 
#' @name dataMultiv
#' @docType data
#' @usage data(dataMultiv)
#' @format This data frame contains the following columns: \describe{
#' \item{PATIENT}{identification of patient} \item{obs}{number of observation
#' for a patient} \item{TIME0}{start of interval} \item{TIME1}{end of interval
#' (death or censoring time)} \item{INDICREC}{recurrent of type 1 status (0:no,
#' 1:yes)} \item{INDICMETA}{recurrent of type 2 status (0:no, 1:yes)}
#' \item{INDICDEATH}{censoring status (0:alive, 1:death)} \item{v1}{dichotomous
#' covariate (0,1)} \item{v2}{dichotomous covariate (0,1)}
#' \item{v3}{dichotomous covariate (0,1)} \item{TIMEGAP}{time to event} }
#' @keywords datasets
NULL





#' Simulated data for recurrent events and a terminal event with weigths using
#' nested case-control design
#' 
#' This contains a simulated sample of of 819 subjects and 1510 observations.
#' This dataset can be used to illustrate how to fit a joint frailty model for
#' data from nested case-control studies.
#' 
#' 
#' @name dataNCC
#' @docType data
#' @usage data(dataNCC)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of patient} \item{cov1}{dichotomous covariate
#' (0,1)} \item{cov2}{dichotomous covariate (0,1)} \item{t.start}{start of
#' interval} \item{t.stop}{end of interval (death or censoring time)}
#' \item{gaptime}{time to event} \item{event}{recurrent event status (0:no,
#' 1:yes)} \item{deathdays}{time of terminal event (death or right-censoring)}
#' \item{death}{censoring status (0:alive, 1:death)} \item{ncc.wts}{weights for
#' NCC design} }
#' @keywords datasets
NULL


#' Simulated data with two levels of grouping
#' 
#' This contains a simulated sample of 400 observations which allow
#' establishing 20 clusters with 4 subgroups and 5 subjects in each subgroup,
#' in order to obtain two levels of grouping. This data set is useful to
#' illustrate how to fit a nested model. Two independent gamma frailty
#' parameters with a variance fixed at 0.1 for the cluster effect and at 0.5
#' for the subcluster effect were generated. Independent survival times were
#' generated from a simple Weibull baseline risk function. The percentage of
#' censoring data was around 30 per cent. The right-censoring variables were
#' generated from a uniform distribution on [1,36] and a left-truncating
#' variable was generated with a uniform distribution on [0,10]. Observations
#' were included only if the survival time is greater than the truncated time.
#' 
#' 
#' @name dataNested
#' @docType data
#' @usage data(dataNested)
#' @format This data frame contains the following columns: \describe{
#' \item{group}{group identification variable} \item{subgroup}{subgroup
#' identification variable} \item{t1}{start of interval (0 or truncated time)}
#' \item{t2}{end of interval (death or censoring time)} \item{event}{censoring
#' status (0: alive, 1: death)} \item{cov1}{dichotomous covariate (0,1)}
#' \item{cov2}{dichotomous covariate (0,1)} }
#' @source
#' 
#' V. Rondeau, L. Filleul, P. Joly (2006). Nested frailty models using maximum
#' penalized likelihood estimation. \emph{Statistics in Medecine}, \bold{25},
#' 4036-4052.
#' @keywords datasets
NULL



#' Rehospitalization colorectal cancer
#' 
#' This contains rehospitalization times after surgery in patients diagnosed
#' with colorectal cancer
#' 
#' 
#' @name readmission
#' @docType data
#' @usage data(readmission)
#' @format This data frame contains the following columns: \describe{
#' \item{id}{identification of each subject. Repeated for each recurrence}
#' \item{enum}{which readmission} \item{t.start}{start of interval (0 or
#' previous recurrence time)} \item{t.stop}{recurrence or censoring time}
#' \item{time}{interocurrence or censoring time} \item{event}{rehospitalization
#' status. All event are 1 for each subject excepting last one that it is 0}
#' \item{chemo}{Did patient receive chemotherapy? 1: No; 2:Yes}
#' \item{sex}{gender: 1:Males 2:Females} \item{dukes}{Dukes' tumoral stage:
#' 1:A-B; 2:C 3:D} \item{charlson}{Comorbidity Charlson's index. Time-dependent
#' covariate.  0: Index 0; 1: Index 1-2; 3: Index >=3 } \item{death}{death
#' indicator. 1:dead and 0:alive } }
#' @source
#' 
#' Gonzalez, JR., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M.,
#' Cambray, M. and Borras, JM (2005). Sex differences in hospital readmission
#' among colorectal cancer patients. \emph{Journal of Epidemiology and
#' Community Health}, \bold{59}, 6, 506-511.
#' @keywords datasets
NULL





#' Transformed Readmission Data for Illness-Death Modeling
#'
#' A dataset derived from the \code{readmission} data (originally from the
#' \code{frailtypack} package, related to rehospitalization times after surgery
#' in colorectal cancer patients). This transformed version reshapes the data
#' to fit a standard illness-death model framework, focusing on the first
#' event (rehospitalization) and the terminal event (death). Recurrent
#' rehospitalization events beyond the first one are excluded. Time is scaled to years.
#'
#' @name readmission2
#' @docType data
#' @usage data(readmission2)
#' @format A data frame with one row per subject, containing columns suitable
#'   for use with the \code{IllnessDeath} function:
#'   \describe{
#'     \item{id}{Unique subject identification number.}
#'     \item{observed_disease_time}{Time (in years since surgery) to either the first rehospitalization (illness), death, or administrative censoring, whichever occurred first.}
#'     \item{observed_death_time}{Time (in years since surgery) to either death or administrative censoring.}
#'     \item{disease_status}{Indicator for the non-terminal event (first rehospitalization). 1 if the subject experienced a first rehospitalization before death/censoring, 0 otherwise.}
#'     \item{death_status}{Indicator for the terminal event (death). 1 if the subject died, 0 if censored.}
#'     \item{dukes}{Dukes' tumoral stage at baseline (Factor or numeric: 1:A-B; 2:C; 3:D).}
#'     \item{sex}{Gender (Factor or numeric: 1:Male; 2:Female).}
#'     \item{charlson}{Comorbidity Charlson's index at baseline (Factor or numeric: 0: Index 0; 1: Index 1-2; 3: Index >=3). Note: Original data had this as time-dependent, this version likely uses the baseline value.}
#'     \item{chemo}{Indicator whether patient received chemotherapy (Factor or numeric: 1:No; 2:Yes).}
#'     \item{group}{An example grouping variable (numeric, derived from id mod 10 + 1), useful for fitting grouped frailty models.}
#' }
#' @details
#' The transformation process involved:
#' \enumerate{
#'   \item Starting with the original \code{readmission} data.
#'   \item Excluding recurrent rehospitalization events, keeping only the interval from surgery (t.start=0) to the first event (event=1) or censoring (event=0).
#'   \item Reshaping the data so each row represents one subject.
#'   \item Defining \code{observed_disease_time} and \code{disease_status} based on the first event interval (t.stop when t.start=0).
#'   \item Defining \code{observed_death_time} and \code{death_status} based on the overall follow-up time and final death status for the subject. If a subject had a first event and then further follow-up, the death time comes from the second interval if available.
#'   \item Scaling time variables (\code{t.stop}) from days (assumed) to years by dividing by 365.
#'   \item Copying baseline covariates (dukes, sex, charlson, chemo) from the subject's first record.
#' }
#' This dataset is intended primarily for demonstrating the \code{IllnessDeath} function.
#'
#' @source Derived from the \code{readmission} dataset, originally described in:
#' Gonzalez, JR., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M.,
#' Cambray, M. and Borras, JM (2005). Sex differences in hospital readmission
#' among colorectal cancer patients. \emph{Journal of Epidemiology and
#' Community Health}, \bold{59}, 6, 506-511.
#' @keywords datasets
NULL




#' Filtered Paquid Sample Data Set for Illness-Death Modeling
#'
#' A dataset derived from \code{Paq1000} (originally a sample from the Paquid
#' study available via the \code{SmoothHazard} package). This version excludes
#' subjects where the age of dementia diagnosis or censoring ('r') was exactly 
#' equal to the age at study entry ('e'), ensuring valid observation intervals
#' when using left truncation.
#'
#' @name Paq810
#' @docType data
#' @usage data(Paq810) 
#'
#' @format A data frame with approximately 810 rows (original 1000 minus exclusions)
#'   and the following 8 columns:
#'   \describe{
#'     \item{dementia}{Dementia status indicator: 0 = non-demented, 1 = demented.}
#'     \item{death}{Death status indicator: 0 = alive, 1 = dead.}
#'     \item{e}{Age at entry into the study (left-truncation time).}
#'     \item{l}{For demented subjects: age at the visit *before* the diagnostic visit. For non-demented subjects: age at the last visit.}
#'     \item{r}{Age at dementia diagnosis or censoring for dementia (event/censoring time for 0->1 transition). Guaranteed to be > `e`.}
#'     \item{t}{Overall exit age. For dead subjects: age at death. For alive subjects: age at the latest news (censoring time for death).}
#'     \item{certif}{Primary school certificate indicator: 0 = without certificate, 1 = with certificate.}
#'     \item{gender}{Gender indicator: 0 = female, 1 = male.}
#' }
#'
#' @details
#' This dataset was created by filtering the \code{Paq1000} data:
#' \code{Paq810 <- Paq1000[Paq1000$r > Paq1000$e, ]}.
#' This step is necessary to prevent issues with \code{survival::Surv(e, r, dementia)}
#' which requires the stop time ('r') to be strictly greater than the start/entry
#' time ('e').
#'
#' The time variables `e`, `l`, `r`, and `t` are all ages in years.
#' This dataset is suitable for fitting illness-death models with left truncation
#' using functions like \code{IllnessDeath}.
#'
#' @source Derived from the \code{Paq1000} dataset, which originates from the
#'   Paquid study and is included in the \code{SmoothHazard} package.
#' @seealso \code{Paq1000}, The \code{SmoothHazard} package.
#' @keywords datasets
NULL




#' Transformed Bone Marrow Transplant Data for Competing Risks
#'
#' A dataset derived from the \code{bmtcrr} data (originally from the
#' \code{casebase} package). This version adds unique subject IDs, an example
#' grouping variable, and calculates an `observed_time` based on age at transplant
#' plus follow-up time in years, setting age 0 as the origin. It's prepared
#' for competing risks analyses, potentially with frailties or left truncation.
#'
#' @name CPRSKbmtcrr
#' @docType data
#' @usage data(CPRSKbmtcrr) 
#'
#' @format A data frame with 177 observations and the following columns:
#'   \describe{
#'     \item{id}{Unique subject identification number.}
#'     \item{Sex}{Gender of the individual (Factor: Male, Female).}
#'     \item{D}{Disease type: ALL or AML (Factor: ALL, AML).}
#'     \item{Phase}{Phase at transplant (Factor: Relapse, CR1, CR2, CR3).}
#'     \item{Age}{Age in years at transplant (start of follow-up).}
#'     \item{Status}{Status indicator: 0 = censored, 1 = relapse, 2 = competing event.}
#'     \item{Source}{Source of stem cells (Factor: BM+PB, PB).}
#'     \item{ftime}{Original failure time in months since transplant.}
#'     \item{group}{Example grouping variable (numeric, derived from id mod 10 + 1).}
#'     \item{observed_time}{Time in years since birth (Age + ftime/12) representing the
#'      time of event or censoring relative to birth as origin.}
#' }
#'
#' @details
#' This dataset was created by taking the original \code{bmtcrr} data from the
#' \code{casebase} package and applying the following transformations:
#' \enumerate{
#'   \item Added a unique subject identifier \code{id}.
#'   \item Added an example grouping variable \code{group} based on \code{id}.
#'   \item Calculated \code{observed_time = Age + ftime/12} to represent the
#'       subject's age at event or censoring, potentially for use with left
#'       truncation at \code{Age}.
#' }
#' The primary event is typically relapse (Status=1), with death without relapse
#' (Status=2) as a competing event. Censoring is Status=0. Note that the time
#' scale for \code{observed_time} is years since birth.
#' 
#'
#' @source Derived from the \code{bmtcrr} dataset available in the \code{casebase} package.
#' @references Scrucca L, Santucci A, Aversa F. Competing risk analysis using R:
#'   an easy guide for clinicians. \emph{Bone Marrow Transplant}. 2007 Aug;40(4):381-7.
#'   doi:10.1038/sj.bmt.1705727.
#' @seealso \code{bmtcrr}. The \code{casebase} package.
#' @keywords datasets
NULL




##' Advanced Ovarian Cancer dataset
##' 
##' This dataset combines the data  that were collected in four double-blind randomized
##' clinical trials in advanced ovarian cancer. In these trials, the objective was to 
##' examine the efficacy of cyclophosphamide plus cisplatin (CP) versus cyclophosphamide 
##' plus adriamycin plus cisplatin (CAP) to treat advanced ovarian cancer. The candidate 
##' surrogate endpoint \bold{S} is progression-free survival time, defined as the time (in years)
##' from randomization to clinical progression of the disease or death. The true endpoint
##' \bold{T} is survival time, defined as the time (in years) from randomization to death of any 
##' cause
##' 
##' 
##' @name dataOvarian
##' @docType data
##' @usage data(dataOvarian)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{patientID}{The identification number of a patient} 
##' \item{trialID}{The center in which a patient was treated}
##' \item{trt}{The treatment indicator, coded as 0 = cyclophosphamide plus cisplatin (CP)
##' and 1 = cyclophosphamide plus adriamycin plus cisplatin(CAP)}
##' \item{timeS}{The candidate surrogate (progression-free survival)}
##' \item{statusS}{Censoring indicator for for Progression-free survival}
##' \item{timeT}{The true endpoint (survival time)}
##' \item{statusT}{Censoring indicator for survival time}
##' }
##' @source
##' 
##' Ovarian cancer Meta-Analysis Project (1991). Cyclophosphamide plus cisplatin plus adriamycin
##' versus Cyclophosphamide, doxorubicin, and cisplatin chemotherapy of ovarian carcinoma: 
##' A meta-analysis. \emph{Classic Papers and Current Comments}, \bold{3}, 237-234.
##' @keywords datasets
NULL

##' Advanced Gastric Cancer dataset
##' 
##' This meta-analysis was carried out by the GASTRIC (Global Advanced/Adjuvant Stomach 
##' Tumor Research international Collaboration) group, using individual data on patients 
##' with curatively resected gastric cancer. Data from all published randomized trials, 
##' with a patient recruitment end date before 2004, and comparing adjuvant chemotherapy 
##' with surgery alone for resectable gastric cancers, were searched electronically. 
##' The candidate surrogate endpoint \bold{S} was Desease-free survival time, defined 
##' as the time (in days) to relapse, second cancer or dead from any cause. The true 
##' endpoint \bold{T} was the overall survival time, defined as the time (in days) from 
##' randomization to death of any cause or to the last follow-up.
##' 
##' 
##' @name gastadj
##' @docType data
##' @usage data(gastadj)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{trialID}{The trial in which the patient was treated}
##' \item{patientID}{The identification number of a patient} 
##' \item{trt}{The treatment indicator, coded as 0 = Control and 1 = Experimental}
##' \item{timeS}{The candidate surrogate (progression-free survival in days) }
##' \item{statusS}{Censoring indicator for for Progression-free survival 
##' (0 = alive and progression-free, 1 = with progression or dead)}
##' \item{timeT}{The true endpoint (overall survival time in days)}
##' \item{statusT}{Censoring indicator for survival time (0 = alive, 1 = dead)}
##' }
##' @source
##' 
##' Oba K, Paoletti X, Alberts S, Bang YJ, Benedetti J, Bleiberg H, Catalona P, 
##' Lordick F, Michiels S, Morita A, Okashi Y, Pignon JP, Rougier P, Sasako M, 
##' Sakamoto J, Sargent D, Shitara K, Van Cutsem E, Buyse M, Burzykowski T on 
##' behalf of the GASTRIC group (2013). Disease-Free Survival as a Surrogate 
##' for Overall Survival in Adjuvant Trials of Gastric Cancer: A Meta-Analysis. 
##' \emph{JNCI: Journal of the National Cancer Institute};\bold{105(21)}:1600-1607
##' @keywords datasets
NULL

##' Longitudinal semicontinuous biomarker dataset (TPJM)
##' 
##' This is a simulated dataset used to illustrate the two-part 
##' joint model included in the longiPenal function.
##' 
##' @name longDat
##' @docType data
##' @usage data(longDat)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{id}{The identification number of a patient}
##' \item{timej}{The measurement times of the biomarker}
##' \item{trtY}{Treatment covariate}
##' \item{Y}{Biomarker value}
##' }
NULL

##' Survival dataset (TPJM)
##' 
##' This is a simulated dataset used to illustrate the two-part 
##' joint model included in the longiPenal function.
##' 
##' @name survDat
##' @docType data
##' @usage data(survDat)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{id}{The identification number of a patient}
##' \item{deathTimes}{The event times (death or censoring)}
##' \item{d}{Censoring indicator}
##' \item{trt}{Treatment covariate}
##' }
NULL

##' Delirium in critically ill ICU patients dataset: the REDUCE clinical trial
##' 
##' This dataset contains an extract of 500 randomly selected patients from the randomized, 
##' double-blind, placebo-controlled REDUCE trial for critically ill patient admited to ICU. 
##' This trial investigated whether Haloperidol (1 or 2 mg) administered as a prophylactic improved 
##' 28-day survival compared to placebo. Recurrent episodes of delirium are 
##' recorded and patients and patients can be censored by death or discharge from the ICU. 
##' 
##' @name reduce
##' @docType data
##' @usage data(reduce)
##' @format This data frame contains the following columns: 
##' \describe{
##' \item{id}{Identification number of a patient} 
##' \item{t.start}{Start time of the interval (0 or time of last recurrence)}
##' \item{t.stop}{Stop time of the interval, either delirium recurrence time or censoring time.}
##' \item{del}{Delirium status}
##' \item{death}{Death status}
##' \item{discharge}{Discharge status}
##' \item{treatment}{Treatment indicator, 
##' 1 if patient was randomized to receive 2mg of Haloperidol, 
##' 0 for control }
##' }
##' @source
##' Van Den Boogaard, M., Slooter, A. J., Bruggemann, R. J., 
##' Schoonhoven, L., Beishuizen, A., Vermeijden, J. W., et al. (2018). 
##' Effect of haloperidol on survival among critically ill adults with a 
##' high risk of delirium: the REDUCE randomized clinical trial. \emph{Jama}, \bold{319(7)}, 680-690.
##' @keywords datasets
NULL