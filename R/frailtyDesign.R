################################################################################
##                         Common Internal Helpers                            ##
################################################################################

# => find the non-centrality parameter for a chi-square statistic
#' @usage NULL
.chi2ncp <- function(alpha, beta, df = 1) {
  cv <- qchisq(alpha, df = df, lower.tail = FALSE)
  ncpmax <- 1
  while (pchisq(cv, df = df, lower.tail = FALSE, ncp = ncpmax) <
    beta) {
    ncpmax <- ncpmax * 2
  }
  return(uniroot(function(ncp) {
    return(pchisq(cv,
      df = df, lower.tail = FALSE,
      ncp = ncp
    ) - beta)
  }, c(0, ncpmax))$root)
}

# Gauss-Laguerre quadrature constants
.glag_xi <- c(
  48.02608557, 38.53068331, 31.40751917, 25.62389423, 20.77647890,
  16.65440771, 13.13028248, 10.12022857, 7.56591623, 5.42533663,
  3.66762272, 2.26994953, 1.21559541, 0.49269174, 0.09330781
)
.glag_wi <- c(
  1.600595e-20, 1.483027e-16, 1.456515e-13, 3.921897e-11, 4.227430e-09,
  2.226317e-07, 6.459927e-06, 1.116744e-04, 1.212436e-03, 8.563878e-03,
  4.020686e-02, 1.264258e-01, 2.630276e-01, 3.422102e-01, 2.182349e-01
)






#' Sample Size calculation and Power Analysis using Gamma-Frailty Models
#'
#' @name frailtyDesign
#'
#' @description{
#' A collection of functions to calculate statistical power and required
#' sample sizes for survival analysis using frailty models, specifically
#' the Shared Frailty Model (SFM), Nested Frailty Model (NFM), Joint
#' Frailty Model (JFM), and General Joint Frailty Model (GJFM).
#' 
#' For each frailty model type (denoted by `*`, where `*` corresponds to
#' SFM, NFM, JFM or GJFM), the package provides two distinct functions:
#' 
#' * `*.power`: Computes the statistical power under given study settings.
#' * `*.ssize`: Determines the required sample size needed to achieve
#' a specified target power under given study settings.
#' }
#'
#' @usage
#' ######################################################
#' ## 1. SHARED FRAILTY MODEL (SFM)
#' ######################################################
#'
#' # Compute power for a given sample size in a SFM
#' # --------------------------------------------------
#' SFM.power(
#'   Groups = 80, ni = 8, ni.type = "max", Acc.Dur = 0,FUP = 12,
#'   FUP.type = "UpToEnd", median.H0 = 1, beta.H0 = 0, beta.HA = log(0.75),
#'   shape.W = 1, theta = 0.25, ratio = 1, samples.mc = 1e4, seed = 42,
#'   timescale = "gap", data.type = "grouped",
#'   cens.par = 5, cens.type = "Expo", statistic = "Wald",
#'   typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' # Compute sample size for a given power in a SFM
#' # --------------------------------------------------
#' SFM.ssize(
#'   power = 0.8, ni = 8, ni.type = "max", Acc.Dur = 0,FUP = 12,
#'   FUP.type = "UpToEnd", median.H0 = 1, beta.H0 = 0, beta.HA = log(0.75),
#'   shape.W = 1, theta = 0.25, ratio = 1, samples.mc = 1e4, seed = 42,
#'   timescale = "gap", data.type = "grouped",
#'   cens.par = 5, cens.type = "Expo", statistic = "Wald",
#'   typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' ######################################################
#' ## 2. NESTED FRAILTY MODEL (NFM)
#' ######################################################
#'
#' # Compute power for a given sample size in a NFM
#' # --------------------------------------------------
#' NFM.power(
#'   Groups = 80, ni = 8, ni.type = "max", kij = 15, kij.type = "max",
#'   Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1,
#'   beta.H0 = 0, beta.HA = log(0.75), shape.W = 1, theta = 0.25, eta = 0.5,
#'   ratio = 1, samples.mc = 1e4, seed = 42,
#'   timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' # Compute sample size for a given power in a NFM
#' # --------------------------------------------------
#' NFM.ssize(
#'   power = 0.8, ni = 8, ni.type = "max", kij = 15, kij.type = "max",
#'   Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1,
#'   beta.H0 = 0, beta.HA = log(0.75), shape.W = 1, theta = 0.25, eta = 0.5,
#'   ratio = 1, samples.mc = 1e4, seed = 42,
#'   timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' ######################################################
#' ## 3. JOINT FRAILTY MODEL (JFM)
#' ######################################################
#'
#' # Compute power for a given sample size in a JFM
#' # --------------------------------------------------
#' JFM.power(
#'   Npts = 400, ni = 8, ni.type = "max", Acc.Dur = 0, FUP = 12,
#'   FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10, betaTest.type = "joint",
#'   betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
#'   shapeR.W = 1, shapeD.W = 1, theta = 0.25, alpha = 1, ratio = 1,
#'   samples.mc = 1e4, seed = 42, timescale = "gap",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' # Compute sample size for a given power in a JFM
#' # --------------------------------------------------
#' JFM.ssize(
#'   power = 0.8, ni = 8, ni.type = "max", Acc.Dur = 0, FUP = 12,
#'   FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10, betaTest.type = "joint",
#'   betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
#'   shapeR.W = 1, shapeD.W = 1, theta = 0.25, alpha = 1, ratio = 1,
#'   samples.mc = 1e4, seed = 42, timescale = "gap",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' ######################################################
#' ## 4. GENERAL JOINT FRAILTY MODEL (GJFM)
#' ######################################################
#'
#' # Compute power for a given sample size in a GJFM
#' # --------------------------------------------------
#' GJFM.power(
#'   Npts = 400, ni = 8, ni.type = "max", Acc.Dur = 0, FUP = 12,
#'   FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
#'   betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75),
#'   betaD.H0 = 0, betaD.HA = log(0.85), shapeR.W = 1, shapeD.W = 1,
#'   theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4,
#'   seed = 42, timescale = "gap",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#'
#' # Compute sample size for a given power in a GJFM
#' # --------------------------------------------------
#' GJFM.ssize(
#'   power = 0.8, ni = 8, ni.type = "max", Acc.Dur = 0, FUP = 12,
#'   FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
#'   betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75),
#'   betaD.H0 = 0, betaD.HA = log(0.85), shapeR.W = 1, shapeD.W = 1,
#'   theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4,
#'   seed = 42, timescale = "gap",
#'   statistic = "Wald", typeIerror = 0.05, test.type = "2-sided"
#' )
#' @param Groups Only in SFM and NFM: A numeric value, where interpretation
#' depends on the \code{data.type} parameter and on the model:
#'   * For SFM, it corresponds to either the number of groups (grouped data) 
#' or the number of subjects (recurrent events data).
#'   * For NFM, it corresponds to either the number of groups (grouped and 
#' recurrent event data) or the number of subjects (multi-type recurrent 
#' events data).
#' 
#' Default is 80.
#' @param ni A numeric value or an array (dim = 2), representing
#' expected values or distribution parameters. Interpretation depends on
#' the \code{data.type} parameter and on the model:
#'   * For SFM, it corresponds to either the expected number of subjects per group (grouped data)
#' or the expected number of recurrent events per subject (recurrent events data).
#'   * For NFM, it corresponds to the expected number of subgroups within each group (grouped data),
#' the expected number of recurrent events per group (recurrent event data), or the number of
#' distinct recurrent event type (multi-type recurrent event data).
#'   * For JFM/GJFM, it corresponds to the expected number of recurrent events per subject.
#'
#' The default value is 8.
#' @param ni.type Character value, specifying \code{ni}. Valid options:
#'   * \code{"max"}: \code{ni} is a fixed number.
#'   * \code{"pois"}: \code{ni} is a mean (parameter of a Poisson distribution).
#'   * \code{"unif"}: \code{ni} is the lower and upper bound parameters of a uniform distribution.
#'
#' Options \code{"pois"} and \code{"unif"} can only be selected for recurrent event data; see Note for details.
#' Default is \code{"max"}.
#' @param median.H0 Only in SFM and NFM: A positive numeric value, used for the scale parameter (Weibull) calculation.
#' If recurrent event data, it is the median gap time to event under the null (excluding censoring times).
#' If grouped data, it is the median time to an event under the null (excluding censoring times).
#' Default is 1.
#' @param Acc.Dur Non-negative numeric value. Parameter for a uniform accrual
#' from time 0 to time \code{Acc.Dur}. Default is 0.
#' @param FUP A positive numeric value of follow-up duration as defined in the
#' study protocol (i.e. administrative censoring). Default is 12.
#' @param FUP.type Character value, indicating the type of follow-up. Valid options:
#'   * \code{"Fixed"}: each subject is followed exactly for \code{FUP} time units after enrollment.
#'   * \code{"UptoEnd"}: global study cutoff at time \code{FUP}; individual follow-up for at most \code{FUP}.
#'
#' Default is \code{"Fixed"}.
#' @param beta.H0 Only in SFM and NFM: log-hazard ratio parameter under the null
#' hypothesis (H0). Default is 0.
#' @param beta.HA Only in SFM and NFM: log-hazard ratio parameter under the
#' alternative hypothesis (HA). Default is log(0.75).
#' @param shape.W Only in SFM and NFM: A positive numeric value, corresponding to
#' the shape parameter (Weibull) of the baseline hazard of the recurrent event. Default is 1.
#' @param theta A positive numeric value, corresponding to the Gamma-frailty
#' variance for the main random effect. Default is 0.25.
#' @param ratio A positive numeric value, corresponding to the allocation ratio
#' (\emph{experimental : control}). Default is 1.
#' @param samples.mc A positive numeric value, corresponding to the number of
#' Monte Carlo samples used to approximate the Fisher information matrix.
#' Default is 1e4.
#' @param cens.par Only in SFM and NFM: A numeric value corresponding
#' to the parameter of the distribution for non-administrative censoring.
#' Default is 10000.
#' @param cens.type Only in SFM and NFM: Character value, specifying
#' the distribution for non-administrative censoring. Valid options:
#'   * \code{"Expo"}: in this case, \code{cens.par} is the median
#'     from an exponential distribution.
#'   * \code{"Unif"}: in this case, \code{cens.par} is the lower and upper
#'      bound parameters of a uniform distribution.
#'   
#' Default is "Expo".
#' @param statistic Type of test statistic used. Currently, only \code{"Wald"} is available.
#' @param typeIerror A numeric value corresponding to the type I error level. Default is 0.05.
#' @param test.type Character value indicating whether It is a one-tailed or two-tailed test.
#'  Valid options are either \code{"1-sided"} or \code{"2-sided"}. Default is "2-sided".
#' @param timescale Character value indicating the timescale when recurrent event data type
#' is considered. Can be either 'gap' or 'calendar'. See note for more detail. Default is 'gap'.
#' @param data.type Only in SFM and NFM: Character value indicating what kind of data we
#' want to consider for the current frailty model. Valid options differ depending on the model:
#'  * For SFM, can be either \code{"grouped"} (corresponding to subjects included in a group)
#'  or \code{"rec_event"} (corresponding to subjects experiencing recurrent events).
#'  * For NFM, the hierarchical structure of the data can be either \code{"grouped"}
#'  (where subjects are included into subgroup and subgroups into groups),
#'  \code{"rec_event1"} (where the group level corresponds to a group (e.g., hospitals)
#'  and subgroup level to a subject) or \code{"rec_event2"} (where the group
#'  level corresponds to a subject and subgroup level to a type of recurrent event).
#'
#' Default is "grouped".
#' @param seed Integer number for random-number generation seed. Ensures reproducibility
#' of the Monte-Carlo simulations. Default is 42.
#' @param kij Only in NFM: A numeric value or an array (dim = 2), representing
#' expected values or distribution parameters. Interpretation depends on the \code{data.type} parameter:
#'  * For grouped data: It is the number of observations per subgroup.
#'  * For recurrent events data: It is the number of observation per subjects.
#'  * For multi-type recurrent events data: It is the number of recurrences for
#' each distinct type of event.
#'
#' Default is 15.
#' @param kij.type Character value, specifying \code{kij}. Valid options:
#'   * \code{"max"}: \code{kij} is a fixed number.
#'   * \code{"pois"}: \code{kij} is a mean (parameter of a Poisson distribution).
#'   * \code{"unif"}: \code{kij} is the lower and upper bound parameters of a uniform distribution.
#'
#' Options \code{"pois"} and \code{"unif"} can only be selected for recurrent event data; see Note for details.
#' Default is \code{"max"}.
#' @param eta Only in NFM and GJFM: positive numeric value, corresponding to an
#' additional Gamma-frailty variance parameter for second-level nesting (NFM)
#' or inter-recurrence dependence (GJFM). Default is 0.5.
#' @param Npts Only in JFM and GJFM: positive numeric value, corresponding to the
#' total number of subjects. Default is 400.
#' @param medianR.H0 Only in JFM and GJFM: positive numeric value, corresponding
#' to the expected median time between two recurrent events under the null (H0),
#'  for the scale parameter (Weibull) calculation. Default is 3.
#' @param medianD.H0 Only in JFM and GJFM: positive numeric value, corresponding
#' to the expected median time to the terminal event under the null (H0),
#'  for the scale parameter (Weibull) calculation. Default is 10.
#' @param power Numeric in (0,0.99]. The target power \eqn{1 - \beta}. Default is 0.8.
#' @param betaTest.type Only in JFM and GJFM: character value indicating which
#' hypothesis is tested when computing power. Our implementation allows either
#' power calculation or sample-size estimation, testing recurrent events alone,
#' terminal event alone or both. Valid options:
#'  \code{"joint"} (for testing both \eqn{\beta_R} and \eqn{\beta_D}), \code{"betaRtest"}
#'  (for testing only \eqn{\beta_R}) or \code{"betaDtest"}
#'  (for testing only \eqn{\beta_D}). Default is \code{"joint"}.
#' @param betaR.H0 Only in JFM and GJFM: numeric value, corresponding to the
#' log-hazard ratios for recurrent events under the null hypothesis (H0). Default is 0.
#' @param betaR.HA Only in JFM and GJFM: numeric value, corresponding to the
#' log-hazard ratios for recurrent events under the alternative hypothesis (HA). Default is log(0.75).
#' @param betaD.H0 Only in JFM and GJFM: numeric value, corresponding to the
#' log-hazard ratios for terminal events under the null hypothesis (H0). Default is 0.
#' @param betaD.HA Only in JFM and GJFM: numeric value, corresponding to the
#' log-hazard ratios for terminal events under the alternative hypothesis (HA). Default is log(0.85).
#' @param shapeR.W Only in JFM and GJFM: positive numeric value, corresponding to
#'  the shape parameter (Weibull) of the recurrent-event baseline hazard function. Default is 1.
#' @param shapeD.W Only in JFM and GJFM: positive numeric value, corresponding to
#'  the shape parameter (Weibull) of the terminal-event baseline hazard function. Default is 1.
#' @param alpha Only in JFM: numeric value, corresponding to the parameter \eqn{\alpha}
#' that modulates the association between recurrent and terminal events. Default is 1.
#'
#' @details{
#' See Dinart et al. (2024) for the original article. We present here the case
#' where we want to assess the treatment effect only. Our null hypothesis is that
#' there is no treatment effect (i.e. zero log-hazard ratio).
#'
#' This approach relies on the squared Wald test to assess the presence of a
#' treatment effect using an estimator \eqn{\hat{\Theta}}. Specifically, under
#' our null hypothesis, the test statistic
#' \eqn{X_w = Z^2 = (\hat{\Theta} - \Theta)^2 / \mathcal{I}^{-1}(\hat{\Theta})}
#' follows a central \eqn{\chi^2_1} distribution (i.e., non-centrality parameter
#' \eqn{\mu = 0}), whereas under the alternative hypothesis, it follows a
#' non-central \eqn{\chi^2_1(\mu)} distribution with \eqn{\mu > 0}.
#'
#' The parameter \eqn{\mu} is estimated algorithmically, and the Fisher
#' information \eqn{\mathcal{I}_1(\hat{\Theta})} is obtained by simulation,
#' leveraging the law of large numbers. Concretely, for an \eqn{M}-sample
#' generated by simulation, the matrix \eqn{\mathcal{I}(\hat{\Theta})} is
#' approximated via the empirical mean of the products
#' \eqn{\partial_{\Theta_k} l(\Theta(i)) \times \partial_{\Theta_l} l(\Theta(i))^\top},
#' \eqn{i \in [\![1, \dots, M]\!]}. The algorithmic estimation for \eqn{\mu}
#' follows the three-step procedure described by Dinart et al. (2024):
#'
#' 1. Compute the \eqn{\alpha}-quantile (denoted \eqn{q_{1,\alpha}}) of a central
#'    chi-square distribution with 1 degree of freedom (\eqn{\chi^2_1}),
#'    given a specified type I error rate \eqn{\alpha}.
#'
#' 2. Determine a non-centrality parameter \eqn{\vartheta} such that
#'    \eqn{1 - P\bigl(\chi^2_1(\vartheta) < q_{1,\alpha}\bigr) > 1 - \beta},
#'    where \eqn{1 - \beta} represents the desired statistical power and
#'    \eqn{1 - P\bigl(\chi^2_1(\vartheta) < q_{1,\alpha}\bigr)} the computed power.
#'
#' 3. Optimize \eqn{\mu} to find the smallest value satisfying that condition for
#' all
#'    \eqn{x \in [0,\vartheta]}, i.e.,
#'    \deqn{
#'      \mu
#'      = \min_{x \in [0,\vartheta]} \Bigl\{ 1 - P\bigl(\chi^2_1(x) < q_{1,\alpha}\bigr)
#'      - (1 - \beta)\Bigr\}.
#'    }
#'
#' Once \eqn{\mu} is estimated, the sample size \eqn{n} is derived from
#' \deqn{
#'   n \,\geq\, \mu
#'   \times \Bigl(\Theta_A^2 \times \mathcal{I}_1(\hat{\Theta})\Bigr)^{-1},
#' }
#' where \eqn{\Theta_A} denotes the parameter value under \eqn{H_A}.
#' If we are interested in the evaluation of the power, we estimate the
#' non-centrality parameter under \eqn{H_A} for a given sample size \eqn{N}, then compute
#' the power as \eqn{P\bigl(\chi^2_1(\vartheta) > q_{1,\alpha}\bigr\vert H_A)}.
#'
#' For both the joint frailty model and general joint frailty model, by following
#' the same methodology as in the univariate case, we can derive an expression
#' for the sample size from the generalized Wald statistic. Let \eqn{H_0:
#' (\beta_R = 0) \text{ and } (\beta_D = 0)} vs. \eqn{H_A: (\beta_R = \beta_R^A)
#' \text{ or } (\beta_D = \beta_D^A)}, be our null and alternative hypotheses respectively.
#' This multivariate test then follows a \eqn{\chi^2_Q} distribution, where
#' \eqn{Q} is the rank of the matrix \eqn{C}, corresponding to the number of
#' constraints applied on the parameters under the null hypothesis. The test
#' statistic is:
#' \deqn{
#'   X_W \;=\; n\,\bigl(C\,\Omega\bigr)^\top
#'   \Bigl(C\,\mathcal{I}_1^{-1}(\Omega)\,C^\top\Bigr)^{-1}\bigl(C\,\Omega\bigr)
#'   \;\sim\;\chi^2_q(\mu),
#' }
#'
#' where \eqn{\Omega^\top} is the vector parameter from the corresponding model.
#' From this, we derive a sample size formula:
#' \deqn{
#'   n \,\geq\, \mu \,\Bigl(\bigl(C\,\Omega\bigr)^\top
#'   \bigl(C\,\mathcal{I}_1^{-1}(\Omega)\,C^\top\bigr)^{-1}\bigl(C\,\Omega\bigr)\Bigr)^{-1}.
#' }
#' For instance, for the JFM, we have
#' \eqn{\Omega^\top = (\beta_R, \beta_D, r_0(.), h_0(.), \theta, \alpha)^\top)}.
#' If we want to test a treatment effect on both the recurrent and the terminal event
#' (i.e. \eqn{H_0 : \beta_R=0 \text{ and } \beta_D=0}), hence:
#' \eqn{C\times \Omega = \begin{pmatrix}1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \end{pmatrix}
#' \times (\beta_R, \beta_D, r_0(.), h_0(.), \theta, \alpha)^\top}
#' }
#'
#' @return
#' For SFM and NFM, the *.power() function returns a list with:
#' \itemize{
#'   \item \code{estimated.power}: the estimated power.
#'   \item \code{number.events}: an array (dim = 2) containing the number of recurrent events
#'   under the null and alternative hypotheses, respectively.
#' }
#'
#' The *.ssize() function returns a list with:
#' \itemize{
#'   \item \code{Groups}: the total number of subjects (recurrent event data)
#'   or number of groups (clustered data).
#'   \item \code{number.events}: an array (dim = 2) containing the number of events
#'   under the null and alternative hypotheses, respectively.
#' }
#'
#' For JFM and GJFM, the *.power() function returns a list with:
#' \itemize{
#'   \item \code{estimated.power}: the estimated power.
#'   \item \code{events.rec}: an array (dim = 2) containing the number of recurrent events
#'   under the null and alternative hypotheses, respectively.
#'   \item \code{events.D}: an array (dim = 2) containing the number of terminal events
#'   under the null and alternative hypotheses, respectively.
#' }
#'
#' The *.ssize() function returns a list with:
#' \itemize{
#'   \item \code{Npts}: the computed total number of subjects.
#'   \item \code{events.rec}: an array (dim = 2) containing the number of recurrent events
#'   under the null and alternative hypotheses, respectively.
#'   \item \code{events.D}: an array (dim = 2) containing the number of terminal events
#'   under the null and alternative hypotheses, respectively.
#' }
#'
#' All returned lists additionally include several input parameters: \code{target.power} or
#' \code{Groups}/\code{Npts} (depending on the called function), \code{ni}, \code{FUP},
#' \code{FUP.type}, \code{Acc.Dur}, \code{ratio}, \code{data.type}, the test type \code{testType},
#' \code{alpha}, \code{theta}, \code{eta} (for NFM and GJFM) and \code{samplesMC}. For SFM and NFM,
#' we have the corresponding hazard ratios (\code{HR.H0}, \code{HR.HA}) from the given betas, 
#' \code{median.H0}. For JGM and GJFM, we have the corresponding hazard ratios (\code{HR.R0}, 
#' \code{HR.RA}, \code{HR.D0}, \code{HR.DA}) from the given betas, \code{medianR.H0}, \code{medianD.H0},
#'  and the testing structure \code{tested.structure}. Along with that, \code{model} 
#' ("SFM", "NFM", "JFM" or "GJFM"), \code{method} ("power" or "ssize") and \code{timescale}
#' ("gap" or "calendar") are also included.
#'
#' All these parameters are utilized by the S3 methods \code{print.frailtyDesign} and
#' \code{summary.frailtyDesign} for further detailed output.
#'
#' @note
#' Internally, these functions rely on extensive numerical integration using Gaussian-Laguerre
#' quadrature to approximate the Fisher information. As such, computations may become resource-intensive.
#' You may need to adjust the parameter \code{samples.mc} or other integration parameters to enhance
#' computational performance or precision.
#'
#' In this implementation, users must provide both the median time and the shape parameter explicitly;
#' the scale parameter is then computed automatically. Under the Weibull distribution, the median time
#' \eqn{t_{1/2}} relates to the scale and shape parameters via:
#' \eqn{t_{1/2} = \text{scale} \times \log(2)^{1/\text{shape}}}. Consequently, the scale parameter is
#' calculated as: \eqn{\text{scale} = \frac{t_{1/2}}{\log(2)^{1/\text{shape}}}}.
#'
#' For both SFM and NFM, when analyzing grouped data, the arguments \code{ni.type}
#' and \code{kij.type} are restricted to the value \code{"max"} to define an
#' exact sample size. For instance, specifying \code{ni.type = "pois"} would
#' represent a mean number of subgroups per group, thereby precluding
#' the determination of a precise sample size.
#'
#' In NFM, the parameter \code{"rec_event2"} might initially appear difficult to interpret.
#' As clarified by Derek et al. (2024), multitype recurrent events include situations such as
#' transient ischemic attacks classified by anatomical location in cardiovascular studies, or
#' migraines differentiated according to severity in neurological research.
#'
#' In survival analysis involving recurrent events, the interpretation of
#' regression coefficients (\eqn{\beta}) is contingent upon the chosen timescale.
#' This distinction is crucial, as the timescale directly influences the risk
#' assessment and the corresponding interpretation of model parameters.
#' * When employing a gap timescale, the timescale resets after each event,
#' measuring the duration until the next occurrence. Consequently, the regression
#' coefficients represent the modification of the inter-event risk, reflecting
#' how the treatment influence the hazard of experiencing a subsequent event
#' after the previous one. This approach focuses on the conditional risk between
#' events.
#' * In contrast, utilizing a calendar timescale measures the time from a fixed
#' origin, such as study entry, without resetting after each event. Here, the
#' regression coefficients pertain to the modification of the risk since the
#' initiation of the study, indicating how the treatment affect the hazard of
#' experiencing events over the entire follow-up period. This approach
#' focuses on the cumulative risk from the study entry.
#'
#' @author
#' Original code by Dinart Derek. Implementation by Adrien Orué.
#'
#' @seealso \code{\link[frailtypack]{frailtyPenal}}, \code{\link[frailtypack]{print.frailtyDesign}},
#' \code{\link[frailtypack]{summary.frailtyDesign}}
#'
#' @references
#' Derek Dinart, Carine Bellera & Virginie Rondeau (09 Feb 2024).
#'         Sample size estimation for recurrent event data using multifrailty
#'         and multilevel survival models,
#'         \emph{Journal of Biopharmaceutical Statistics},
#'         DOI: 10.1080/10543406.2024.2310306.
#'
#' @examples
#' \donttest{
#' # Example 1 (SFM): a total of 400 patients (1:1 randomization scheme),
#' # with a fixed number of 3 recurrent events per patient. Gamma-frailty
#' # variance of 0.5. Expected hazard ratio of 0.7, time-to-death are uniformly
#' # distributed, with a mean time to death of (3+10)/2=6.5 years. Each subject is
#' # followed-up for a maximum of 6 years, with a median time-to-event of 1.5 years.
#' # Patients are recruited over a 0.5-year period.
#' SFM.power(
#'   Groups = 400, ni = 3, ni.type = "max",
#'   FUP = 6, Acc.Dur = 0.5, median.H0 = 1.5,
#'   beta.HA = log(0.7), theta = 0.5,
#'   cens.par = c(3, 10), cens.type = "Unif",
#'   data.type = "rec_event"
#' ) # power ~ 90%
#'
#'
#' # Example 2 (NFM): same parameters as above, but we now assume that we have
#' # 40 hospitals, 10 subjects per hospital (10 × 40 = 400 subjects in total)
#' # and 3 recurrent events per subject.
#' NFM.power(
#'   Groups = 40, ni = 10, ni.type = "max", kij = 3, kij.type = "max",
#'   FUP = 6, Acc.Dur = 0.5, median.H0 = 1.5,
#'   beta.HA = log(0.7), theta = 0.5,
#'   cens.par = c(3, 10), cens.type = "Unif",
#'   data.type = "rec_event1"
#' ) # power ~ 83%
#'
#'
#' # Example 3 (NFM): we aim to compute the required sample size to achieve
#' # 80% power for detecting a hazard ratio of 0.75 in a neurological study,
#' # where migraine episodes experienced by subjects are classified into three
#' # severity subtypes (mild, moderate, severe). For each subject, we anticipate
#' # a mean number of 2 migraine episodes per severity subtype, with a median
#' # time-to-event of 6 months. The study duration includes a 1-year accrual period
#' # followed by a 5-year total follow-up. All subjects will be followed until
#' # the end of the study.
#' NFM.ssize(
#'   power = 0.80, ni = 3, ni.type = "max", kij = 2, kij.type = "pois",
#'   FUP = 5, Acc.Dur = 1, FUP.type = "uptoend", median.H0 = 0.5,
#'   beta.HA = log(0.75), data.type = "rec_event2"
#' ) # sample size ~ 363 patients
#'
#' # Example 4 (JFM): power estimation, testing a treatment effect on recurrent
#' # events only. We assume a uniformly distributed number of recurrent events,
#' # ranging from 1 to 6 recurrent events per subject. The allocation ratio
#' # experimental:control is 2:1, and the follow-up is 10 weeks.
#' # The expected hazard ratio is 0.70 for recurrent event and 0.90 for the
#' # terminal event. We have chosen 0.5 as the variance of the frailties.
#' JFM.power(
#'   Npts = 400, ni = c(1, 6), ni.type = "unif",
#'   FUP = 10, FUP.type = "fixed", ratio = 2,
#'   betaTest.type = "betaRtest", betaR.HA = log(.70), betaD.HA = log(.90),
#'   theta = .5
#' ) # power ~ 76%
#'
#'
#' # Example 5 (JFM): sample size calculation, to assess the treatment effect on
#' # both recurrent and terminal events. We want to achieve an 80% power.
#' # We anticipate a maximum of 5 recurrent events, over a 6-year period and a
#' # 0.5-year accrual period. We assume that the gamma-frailty variance is 0.5.
#' # For the control group, we expect a 2-year and a 5-year median time-to-event
#' # for recurrent events and terminal events, respectively. We consider a 30%
#' # and 20% risk reduction for recurrent events and terminal event, respectively.
#' JFM.ssize(
#'   power = 0.80, ni = 9,
#'   FUP = 6, Acc.Dur = 1.5, medianR.H0 = 2, medianD.H0 = 5,
#'   betaTest.type = "joint", betaR.HA = log(.70), betaD.HA = log(.80), theta = .5
#' ) # sample size ~ 445 patients / ~ approx 
#'
#'
#' # Example 6: Sample size calculation for GJFM
#' # Same as above, but with two random effects (with two Gamma-frailty variances
#' # theta and eta). To ensure sample size estimation stability, we use 10000
#' # Monte-Carlo samples.
#' GJFM.ssize(
#'   power = 0.80, ni = 5,
#'   FUP = 6, Acc.Dur = 0.5, medianR.H0 = 2, medianD.H0 = 5,
#'   betaR.HA = log(0.70), betaD.HA = log(0.80), theta = 0.5, eta = 0.75,
#'   samples.mc = 1e5
#' ) # sample size ~ 705 patients / ~ approx 4 min.
#'
#'
#' # Example 7:
#' # --------------------------------------------------
#' # Post-hoc power analysis for a Joint Frailty Model
#' # --------------------------------------------------
#' # See original article by Gonzalez et al. (2005)
#' data(readmission)
#' modJFM <- frailtyPenal(
#'   Surv(time, event) ~ cluster(id) + as.factor(chemo) + terminal(death),
#'   formula.terminalEvent = ~ as.factor(chemo), data = readmission,
#'   hazard = "Weibull"
#' )
#'
#' # Test both recurrent and death events
#' # # - Let us assume an underlying Poisson distribution for ni.type. The
#' # # empirical mean of the number of recurrent events per patients is: ni = 1.136476.
#' # # - For the null hypothesis, let us consider betaR.H0 = betaD.H0 = 0. For the
#' # # alternative hypothesis, we use the estimated parameters for betaD.HA and
#' # # betaR.HA.
#' # # - "Patients were actively followed up until June 2002" -> the follow-up
#' # # type is "UpToEnd".
#' # # - "The study took place in the Hospital de Bellvitge, [...] between
#' # # January 1996 and December 1998" -> the accrual time is approximately 3 years
#' # # - We can assume that the study duration is approximately 6 years
#'
#' ni <- 1.136476
#' ni.type <- "Pois"
#' Acc.Dur <- 3 * 365.25 # time unit = days
#' FUP <- 6 * 365.25     # same as above
#' betaR.HA <- as.numeric(modJFM$coef[1]) # else "Named numeric"
#' betaD.HA <- as.numeric(modJFM$coef[2]) # same as above
#' med <- modJFM$scale.weib * log(2)^(1 / modJFM$shape.weib)
#' medianR.H0 <- med[1]
#' medianD.H0 <- med[2]
#' shapeR.W <- modJFM$shape[1]
#' shapeD.W <- modJFM$shape[2]
#' theta <- modJFM$theta
#' alpha <- modJFM$alpha
#' Npts <- length(unique(readmission[, "id"])) # 403 patients
#' nTreated <- length(unique(readmission[readmission$chemo == "Treated", "id"])) #217 treated patients
#' ratio <- nTreated / (Npts - nTreated)
#'
#' JFM.power(
#'   Npts = Npts, ni = ni, ni.type = ni.type,
#'   Acc.Dur = Acc.Dur, FUP = FUP, medianR.H0 = medianR.H0, medianD.H0 = medianD.H0,
#'   betaTest.type = "joint", betaR.HA = betaR.HA, betaD.HA = betaD.HA,
#'   shapeR.W = shapeR.W, shapeD.W = shapeD.W, theta = theta, alpha = alpha,
#'   ratio = ratio
#' ) # power ~ 92%
#'
#' # --------------------------------------------------
#' # Required sample size under the same setting
#' # --------------------------------------------------
#' # Here, let us consider that readmission is a “pilot study” with 403 patients,
#' # from which we estimate parameters. Under this scenario, let us compute the
#' # needed sample size, but to achieve an 80% power.
#'
#' JFM.ssize(
#'   power = 0.80, ni = ni, ni.type = ni.type,
#'   Acc.Dur = Acc.Dur, FUP = FUP, medianR.H0 = medianR.H0, medianD.H0 = medianD.H0,
#'   betaTest.type = "joint", betaR.HA = betaR.HA, betaD.HA = betaD.HA,
#'   shapeR.W = shapeR.W, shapeD.W = shapeD.W, theta = theta, alpha = alpha,
#'   ratio = ratio
#' ) # 289 patients needed under the same settings vs. 403
#' }
#' 
#' @import MASS
#' @importFrom matrixcalc is.positive.definite is.square.matrix
#' @importFrom stats ave qchisq rexp rpois runif uniroot
#' @keywords methods power
#' @md










################################################################################
##                             Called Functions                               ##
################################################################################
#' @rdname frailtyDesign
#' @usage NULL
#' @export
SFM.power <- function(Groups = 80,
                      ni = 8,
                      ni.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      median.H0 = 1,
                      beta.H0 = 0,
                      beta.HA = log(0.75),
                      shape.W = 1,
                      theta = 0.25,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      data.type = "grouped",
                      cens.par = 5,
                      cens.type = "Expo",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Power_SFM(
    Groups = Groups,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    median.H0 = median.H0,
    beta.H0 = beta.H0,
    beta.HA = beta.HA,
    shape.W = shape.W,
    theta = theta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    data.type = data.type,
    cens.par = cens.par,
    cens.type = cens.type,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  # Append additional fields for summary display
  out$model <- "SFM"
  out$method <- "power"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$samplesMC <- samples.mc
  out$ratio <- ratio
  out$data.type <- data.type
  out$FUP.type <- FUP.type
  out$FUP <- FUP
  out$ni <- ni
  out$ni.type <- ni.type
  out$Acc.Dur <- Acc.Dur
  out$Groups <- Groups
  out$theta <- theta

  out$HR.H0 <- exp(beta.H0)
  out$HR.HA <- exp(beta.HA)

  out$median.H0 <- median.H0


  if (!is.null(out$events)) {
    out$number.events <- ceiling(out$events)
    out$events <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}


#' @rdname frailtyDesign
#' @usage NULL
#' @export
SFM.ssize <- function(power = 0.80,
                      ni = 8,
                      ni.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      median.H0 = 1,
                      beta.H0 = 0,
                      beta.HA = log(0.75),
                      shape.W = 1,
                      theta = 0.25,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      data.type = "grouped",
                      cens.par = 5,
                      cens.type = "Expo",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Nsn_SFM(
    power = power,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    median.H0 = median.H0,
    beta.H0 = beta.H0,
    beta.HA = beta.HA,
    shape.W = shape.W,
    theta = theta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    data.type = data.type,
    cens.par = cens.par,
    cens.type = cens.type,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "SFM"
  out$method <- "ssize"
  out$timescale <- timescale
  out$testType <- test.type
  out$samplesMC <- samples.mc
  out$alpha <- typeIerror
  out$data.type <- data.type
  out$ratio <- ratio
  out$FUP.type <- FUP.type
  out$FUP <- FUP
  out$ni <- ni
  out$ni.type <- ni.type
  out$Acc.Dur <- Acc.Dur
  out$theta <- theta
  
  out$estimated.power <- NULL
  out$target.power <- power

  out$HR.H0 <- exp(beta.H0)
  out$HR.HA <- exp(beta.HA)
  out$median.H0 <- median.H0

  if (!is.null(out$events)) {
    out$number.events <- ceiling(out$events)
    out$events <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}


#' @rdname frailtyDesign
#' @usage NULL
#' @export
NFM.power <- function(Groups = 80,
                      ni = 8,
                      ni.type = "max",
                      kij = 15,
                      kij.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      median.H0 = 1,
                      beta.H0 = 0,
                      beta.HA = log(0.75),
                      shape.W = 1,
                      theta = 0.25,
                      eta = 0.5,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      data.type = "grouped",
                      cens.par = 5,
                      cens.type = "Expo",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Power_NFM(
    Groups = Groups,
    ni = ni,
    ni.type = ni.type,
    kij = kij,
    kij.type = kij.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    median.H0 = median.H0,
    beta.H0 = beta.H0,
    beta.HA = beta.HA,
    shape.W = shape.W,
    theta = theta,
    eta = eta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    data.type = data.type,
    cens.par = cens.par,
    cens.type = cens.type,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "NFM"
  out$method <- "power"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$samplesMC <- samples.mc
  out$ratio <- ratio
  out$data.type <- data.type
  out$FUP.type <- FUP.type
  out$FUP <- FUP
  out$ni <- ni
  out$ni.type <- ni.type
  out$kij <- kij
  out$kij.type <- kij.type
  out$Acc.Dur <- Acc.Dur
  out$theta <- theta
  out$eta <- eta

  out$Groups <- Groups

  out$HR.H0 <- exp(beta.H0)
  out$HR.HA <- exp(beta.HA)
  out$median.H0 <- median.H0

  if (!is.null(out$events)) {
    out$number.events <- ceiling(out$events)
    out$events <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}

#' @rdname frailtyDesign
#' @usage NULL
#' @export
NFM.ssize <- function(power = 0.80,
                      ni = 8,
                      ni.type = "max",
                      kij = 15,
                      kij.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      median.H0 = 1,
                      beta.H0 = 0,
                      beta.HA = log(0.75),
                      shape.W = 1,
                      theta = 0.25,
                      eta = 0.5,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      data.type = "grouped",
                      cens.par = 5,
                      cens.type = "Expo",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Nsn_NFM(
    power = power,
    ni = ni,
    ni.type = ni.type,
    kij = kij,
    kij.type = kij.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    median.H0 = median.H0,
    beta.H0 = beta.H0,
    beta.HA = beta.HA,
    shape.W = shape.W,
    theta = theta,
    eta = eta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    data.type = data.type,
    cens.par = cens.par,
    cens.type = cens.type,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "NFM"
  out$method <- "ssize"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$data.type <- data.type
  out$samplesMC <- samples.mc
  out$ratio <- ratio
  out$FUP.type <- FUP.type
  out$ni <- ni
  out$ni.type <- ni.type
  out$kij <- kij
  out$kij.type <- kij.type
  out$theta <- theta
  out$eta <- eta

  out$FUP <- FUP
  out$Acc.Dur <- Acc.Dur

  out$target.power <- power
  out$estimated.power <- NULL

  out$HR.H0 <- exp(beta.H0)
  out$HR.HA <- exp(beta.HA)
  out$median.H0 <- median.H0

  if (!is.null(out$events)) {
    out$number.events <- ceiling(out$events)
    out$events <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}

#' @rdname frailtyDesign
#' @usage NULL
#' @export
JFM.power <- function(Npts = 400,
                      ni = 8,
                      ni.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      medianR.H0 = 3,
                      medianD.H0 = 10,
                      betaTest.type = "joint",
                      betaR.H0 = 0,
                      betaR.HA = log(0.75),
                      betaD.H0 = 0,
                      betaD.HA = log(0.85),
                      shapeR.W = 1,
                      shapeD.W = 1,
                      theta = 0.25,
                      alpha = 1,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Power_JFM(
    Npts = Npts,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    medianR.H0 = medianR.H0,
    medianD.H0 = medianD.H0,
    betaTest.type = betaTest.type,
    betaR.H0 = betaR.H0,
    betaR.HA = betaR.HA,
    betaD.H0 = betaD.H0,
    betaD.HA = betaD.HA,
    shapeR.W = shapeR.W,
    shapeD.W = shapeD.W,
    theta = theta,
    alpha = alpha,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "JFM"
  out$method <- "power"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$samplesMC <- samples.mc
  out$Npts <- Npts
  out$ratio <- ratio
  out$ni <- ni
  out$ni.type <- ni.type
  out$theta <- theta
 
  out$FUP <- FUP
  out$Acc.Dur <- Acc.Dur

  out$HR.R0 <- exp(betaR.H0)
  out$HR.RA <- exp(betaR.HA)
  out$HR.D0 <- exp(betaD.H0)
  out$HR.DA <- exp(betaD.HA)

  out$medianR.H0 <- medianR.H0
  out$medianD.H0 <- medianD.H0

  out$tested.structure <- betaTest.type

  if (!is.null(out$events.rec)) {
    out$number.events.rec <- ceiling(out$events.rec)
    out$events.rec <- NULL
  }
  if (!is.null(out$events.D)) {
    out$number.events.D <- ceiling(out$events.D)
    out$events.D <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}


#' @rdname frailtyDesign
#' @usage NULL
#' @export
JFM.ssize <- function(power = 0.80,
                      ni = 8,
                      ni.type = "max",
                      Acc.Dur = 0,
                      FUP = 12,
                      FUP.type = "UpToEnd",
                      medianR.H0 = 3,
                      medianD.H0 = 10,
                      betaTest.type = "joint",
                      betaR.H0 = 0,
                      betaR.HA = log(0.75),
                      betaD.H0 = 0,
                      betaD.HA = log(0.85),
                      shapeR.W = 1,
                      shapeD.W = 1,
                      theta = 0.25,
                      alpha = 1,
                      ratio = 1,
                      samples.mc = 1e4,
                      seed = 42,
                      timescale = "gap",
                      statistic = "Wald",
                      typeIerror = 0.05,
                      test.type = "2-sided") {
  out <- Nsn_JFM(
    power = power,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    medianR.H0 = medianR.H0,
    medianD.H0 = medianD.H0,
    betaTest.type = betaTest.type,
    betaR.H0 = betaR.H0,
    betaR.HA = betaR.HA,
    betaD.H0 = betaD.H0,
    betaD.HA = betaD.HA,
    shapeR.W = shapeR.W,
    shapeD.W = shapeD.W,
    theta = theta,
    alpha = alpha,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "JFM"
  out$method <- "ssize"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$samplesMC <- samples.mc
  out$ratio <- ratio
  out$ni <- ni
  out$ni.type <- ni.type  
  out$theta <- theta

  out$target.power <- power
  out$estimated.power <- NULL

  out$HR.R0 <- exp(betaR.H0)
  out$HR.RA <- exp(betaR.HA)
  out$HR.D0 <- exp(betaD.H0)
  out$HR.DA <- exp(betaD.HA)

  out$medianR.H0 <- medianR.H0
  out$medianD.H0 <- medianD.H0

  out$tested.structure <- betaTest.type

  if (!is.null(out$events.rec)) {
    out$number.events.rec <- ceiling(out$events.rec)
    out$events.rec <- NULL
  }
  if (!is.null(out$events.D)) {
    out$number.events.D <- ceiling(out$events.D)
    out$events.D <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}


#' @rdname frailtyDesign
#' @usage NULL
#' @export
GJFM.power <- function(Npts = 400,
                       ni = 8,
                       ni.type = "max",
                       Acc.Dur = 0,
                       FUP = 12,
                       FUP.type = "UpToEnd",
                       medianR.H0 = 3,
                       medianD.H0 = 10,
                       betaTest.type = "joint",
                       betaR.H0 = 0,
                       betaR.HA = log(0.75),
                       betaD.H0 = 0,
                       betaD.HA = log(0.85),
                       shapeR.W = 1,
                       shapeD.W = 1,
                       theta = 0.25,
                       eta = 0.5,
                       ratio = 1,
                       samples.mc = 1e4,
                       seed = 42,
                       timescale = "gap",
                       statistic = "Wald",
                       typeIerror = 0.05,
                       test.type = "2-sided") {
  out <- Power_GJFM(
    Npts = Npts,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    medianR.H0 = medianR.H0,
    medianD.H0 = medianD.H0,
    betaTest.type = betaTest.type,
    betaR.H0 = betaR.H0,
    betaR.HA = betaR.HA,
    betaD.H0 = betaD.H0,
    betaD.HA = betaD.HA,
    shapeR.W = shapeR.W,
    shapeD.W = shapeD.W,
    theta = theta,
    eta = eta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "GJFM"
  out$method <- "power"
  out$timescale <- timescale
  out$testType <- test.type
  out$alpha <- typeIerror
  out$samplesMC <- samples.mc
  out$Npts <- Npts
  out$ratio <- ratio
  out$ni <- ni
  out$ni.type <- ni.type
  out$theta <- theta
  out$eta <- eta

  out$FUP <- FUP
  out$Acc.Dur <- Acc.Dur

  out$HR.R0 <- exp(betaR.H0)
  out$HR.RA <- exp(betaR.HA)
  out$HR.D0 <- exp(betaD.H0)
  out$HR.DA <- exp(betaD.HA)

  out$medianR.H0 <- medianR.H0
  out$medianD.H0 <- medianD.H0

  out$tested.structure <- betaTest.type

  if (!is.null(out$events.rec)) {
    out$number.events.rec <- ceiling(out$events.rec)
    out$events.rec <- NULL
  }
  if (!is.null(out$events.D)) {
    out$number.events.D <- ceiling(out$events.D)
    out$events.D <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}

#' @rdname frailtyDesign
#' @usage NULL
#' @export
GJFM.ssize <- function(power = 0.80,
                       ni = 8,
                       ni.type = "max",
                       Acc.Dur = 0,
                       FUP = 12,
                       FUP.type = "UpToEnd",
                       medianR.H0 = 3,
                       medianD.H0 = 10,
                       betaTest.type = "joint",
                       betaR.H0 = 0,
                       betaR.HA = log(0.75),
                       betaD.H0 = 0,
                       betaD.HA = log(0.85),
                       shapeR.W = 1,
                       shapeD.W = 1,
                       theta = 0.25,
                       eta = 0.5,
                       ratio = 1,
                       samples.mc = 1e4,
                       seed = 42,
                       timescale = "gap",
                       statistic = "Wald",
                       typeIerror = 0.05,
                       test.type = "2-sided") {
  out <- Nsn_GJFM(
    power = power,
    ni = ni,
    ni.type = ni.type,
    Acc.Dur = Acc.Dur,
    FUP = FUP,
    FUP.type = FUP.type,
    medianR.H0 = medianR.H0,
    medianD.H0 = medianD.H0,
    betaTest.type = betaTest.type,
    betaR.H0 = betaR.H0,
    betaR.HA = betaR.HA,
    betaD.H0 = betaD.H0,
    betaD.HA = betaD.HA,
    shapeR.W = shapeR.W,
    shapeD.W = shapeD.W,
    theta = theta,
    eta = eta,
    ratio = ratio,
    samples.mc = samples.mc,
    seed = seed,
    timescale = timescale,
    statistic = statistic,
    typeIerror = typeIerror,
    test.type = test.type
  )

  out$model <- "GJFM"
  out$method <- "ssize"
  out$timescale <- timescale
  out$testType <- test.type
  out$samplesMC <- samples.mc
  out$alpha <- typeIerror
  out$ratio <- ratio
  out$ni <- ni
  out$ni.type <- ni.type
  out$theta <- theta
  out$eta <- eta

  out$FUP <- FUP
  out$Acc.Dur <- Acc.Dur

  out$estimated.power <- NULL
  out$target.power <- power

  out$HR.R0 <- exp(betaR.H0)
  out$HR.RA <- exp(betaR.HA)
  out$HR.D0 <- exp(betaD.H0)
  out$HR.DA <- exp(betaD.HA)

  out$medianR.H0 <- medianR.H0
  out$medianD.H0 <- medianD.H0

  out$tested.structure <- betaTest.type

  if (!is.null(out$events.rec)) {
    out$number.events.rec <- ceiling(out$events.rec)
    out$events.rec <- NULL
  }
  if (!is.null(out$events.D)) {
    out$number.events.D <- ceiling(out$events.D)
    out$events.D <- NULL
  }

  class(out) <- "frailtyDesign"
  return(out)
}



























################################################################################
##                        Computation Functions                               ##
################################################################################

Power_SFM <- function(
    Groups = 400, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1,
    beta.H0 = 0, beta.HA = log(0.75),
    shape.W = 1, theta = 0.25, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo",
    statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(Groups)) {
    stop("Groups is missing")
  }
  if (Groups <= 0) {
    stop("Groups must be strictly positive")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'unif' or 'pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is fixed ")
  }
  if ((!tolower(data.type) %in% c("rec_event", "grouped"))) {
    stop("data.type must be either 'grouped' or 'rec_event' ")
  }
  if (data.type == "grouped") {
    if ((length(ni) == 1) & (ni[1] <= 0)) {
      stop("ni must be strictly positive")
    }
    if ((length(ni) == 2) & ((ni[1] <= 0) | (ni[2] <= 0))) {
      stop("ni must be strictly positive")
    }
  } else {
    if ((length(ni) == 1) & (ni[1] < 0)) {
      stop("ni must be  positive")
    }
    if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
      stop("ni must be  positive")
    }
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be either 'gap' or 'calendar'")
  }
  if (tolower(data.type) == "grouped" & tolower(timescale) != "gap") {
    stop("For grouped data structure, one must use a gap timescale")
  }
  if (tolower(data.type) == "grouped" & tolower(ni.type) != "max") {
    stop("For grouped data structure, ni.type must be max (to ensure a fixed number of subjects per group)")
  }
  if ((!tolower(cens.type) %in% c ("unif", "uniform", "u", "exp", "expo", "exponential"))) {
    stop("cens.type must be either 'Uniform' or 'Exponential ")
  }
  if ((tolower(cens.type) %in% c("u", "unif", "uniform")) & length(cens.par) != 2) {
    stop("cens.par must be a vector of length two if cens.type is uniform ")
  }
  if ((tolower(cens.type) %in% c("expo", "exp", "exponential")) & length(cens.par) != 1) {
    stop("cens.par must be a scalar if cens.type is exponential ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }

  Gc <- Groups
  P_E <- 1 / (1 / ratio + 1)

  theta <- theta

  beta.HA <- beta.HA
  beta.H0 <- beta.H0

  npar <- 4
  if (test.type == "1-sided") {
    alpha <- typeIerror
  }
  if (test.type == "2-sided") {
    alpha <- typeIerror / 2
  }

  scale_W <- median.H0 / log(2)^(1 / shape.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }

  cumsumrev <- function(x) {
    output <- rep(NA, length(x))
    ll <- length(x)
    if (ll > 1) {
      x1 <- x[1:(ll - 1)]
      x2 <- x[2:ll]
      output <- c(x[1], x2 - x1)
    }
    if (ll == 1) {
      output <- x
    }
    return(output)
  }

  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_ni <- 0
  Di_H0 <- Di_HA <- cens_E <- cens_E_ <- 0


  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "fixed"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }
    if (ni_ > 0) {
      if (data.type == "rec_event") {
        Xij <- rep(sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E)), ni_)
        Accrual <- rep(runif(n = 1, min = 0, max = Acc.Dur), ni_)
      } else {
        Xij <- sample(x = c(0, 1), size = ni_, prob = c(1 - P_E, P_E), replace = TRUE)
        Accrual <- runif(n = ni_, min = 0, max = Acc.Dur)
      }

      U <- runif(n = ni_, 0, 1)
      Vi <- rgamma(1, shape = theta^(-1), scale = theta)

      Tij_H0 <- H0_1(s = -log(U) * exp(-(beta.H0 * Xij)) / Vi, y = scale_W, p = shape.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(beta.HA * Xij)) / Vi, y = scale_W, p = shape.W)

      # U1 <-runif(n =1,0,1) ##
      # H0_1(s = -log(U1)*exp(-(beta.HA*Xij))/Vi, y=scale_W,p=shape.W)


      if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
        if (data.type == "rec_event") {
          deathtime <- runif(1, cens.par[1], cens.par[2])
        } else {
          deathtime <- runif(ni_, cens.par[1], cens.par[2])
        }
      } else {
        death_rate <- log(2) / cens.par
        if (data.type == "rec_event") {
          deathtime <- rexp(1, death_rate)
        } else {
          deathtime <- rexp(ni_, death_rate)
        }
      }

      if (tolower(FUP.type) %in% c("f", "fixed")) {
        cens_time <- pmin(FUP, deathtime)
      }
      if (tolower(FUP.type) %in%  c("uptoend", "end")) {
        cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
      }
    } else {
      Xij <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
      Accrual <- runif(n = 1, min = 0, max = Acc.Dur)

      if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
        deathtime <- runif(1, cens.par[1], cens.par[2])
      } else {
        death_rate <- log(2) / cens.par
        deathtime <- rexp(1, death_rate)
      }
      if (FUP.type %in% c("Fixed", "fixed", "f")) {
        cens_time <- pmin(FUP, deathtime)
      }
      if (tolower(FUP.type) %in%  c("uptoend", "end")) {
        cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
      }
      Tij_H0 <- cens_time
      Tij_HA <- cens_time
    }


    if (data.type == "rec_event") {
      Yij_cal_H0 <- pmin(cumsum(Tij_H0), cens_time)
      Yij_cal_HA <- pmin(cumsum(Tij_HA), cens_time)

      if (tolower(timescale) == "gap") {
        Yij_H0 <- cumsumrev(Yij_cal_H0)
        Yij_HA <- cumsumrev(Yij_cal_HA)
      } else {
        Yij_H0 <- Yij_cal_H0
        Yij_HA <- Yij_cal_HA
      }
      Di <- sum(as.numeric((Yij_cal_H0 < cens_time)))
      Di_ <- sum(as.numeric((Yij_cal_HA < cens_time)))
      DiXij <- sum(as.numeric((Yij_cal_H0 < cens_time)) * Xij)
      DiXij_ <- sum(as.numeric((Yij_cal_HA < cens_time)) * Xij)
    } else {
      Yij_H0 <- pmin(Tij_H0, cens_time)
      Yij_HA <- pmin(Tij_HA, cens_time)
      Di <- sum(as.numeric((Yij_H0 < cens_time)))
      Di_ <- sum(as.numeric((Yij_HA < cens_time)))
      DiXij <- sum(as.numeric((Yij_H0 < cens_time)) * Xij)
      DiXij_ <- sum(as.numeric((Yij_HA < cens_time)) * Xij)
    }


    H0Yij_H0 <- sum(c(H0(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))))
    dH0_epsiYij_H0 <- sum(c(dH0_epsi(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))))
    dH0_nuYij_H0 <- sum(c(dH0_nu(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))), na.rm = TRUE)

    H0Yij_HA <- sum(c(H0(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))))
    dH0_epsiYij_HA <- sum(c(dH0_epsi(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))))
    dH0_nuYij_HA <- sum(c(dH0_nu(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))), na.rm = TRUE)


    ln <- c(1:(Di))
    ln_ <- c(1:(Di_))

    score_beta_H0 <- DiXij - (1 / theta + Di) * (theta * sum(H0(t = Yij_H0, y = scale_W, p = shape.W) * Xij *
      exp(beta.H0 * Xij))) / (1 + theta * (H0Yij_H0))


    score_theta_H0 <- log(1 + theta * H0Yij_H0) / theta^2 - (1 / theta + Di) * H0Yij_H0 / (1 + theta * H0Yij_H0) +
      ifelse(Di > 0, sum((Di - ln) / (1 + theta * (Di - ln))), 0)

    score_epsi_H0 <- ifelse((tolower(data.type) == "rec_event"),
      -(1 / theta + Di) * theta * dH0_epsiYij_H0 /
        (1 + theta * H0Yij_H0) - sum(as.numeric((Yij_cal_H0 < cens_time)) * shape.W / scale_W),
      -(1 / theta + Di) * theta * dH0_epsiYij_H0 /
        (1 + theta * H0Yij_H0) - sum(as.numeric((Yij_H0 < cens_time)) * shape.W / scale_W)
    )

    nu_fact_H0 <- ifelse((tolower(data.type) == "rec_event"),
      sum(c(as.numeric((Yij_cal_H0 < cens_time)) * (dh0_nu(t = Yij_H0, y = scale_W, p = shape.W) / h0(t = Yij_H0, y = scale_W, p = shape.W))), na.rm = TRUE),
      sum(as.numeric((Yij_H0 < cens_time)) * (dh0_nu(t = Yij_H0, y = scale_W, p = shape.W) / h0(t = Yij_H0, y = scale_W, p = shape.W)))
    )
    score_nu_H0 <- -(1 / theta + Di) * theta * dH0_nuYij_H0 / (1 + theta * H0Yij_H0) + nu_fact_H0


    score_beta_HA <- DiXij_ - (1 / theta + Di_) * (theta * sum(H0(t = Yij_HA, y = scale_W, p = shape.W) * Xij *
      exp(beta.HA * Xij))) / (1 + theta * (H0Yij_HA))

    score_theta_HA <- log(1 + theta * H0Yij_HA) / theta^2 - (1 / theta + Di_) * H0Yij_HA / (1 + theta * H0Yij_HA) +
      ifelse(Di_ > 0, sum((Di_ - ln_) / (1 + theta * (Di_ - ln_))), 0)

    score_epsi_HA <- ifelse((tolower(data.type) == "rec_event"),
      -(1 / theta + Di_) * theta * dH0_epsiYij_HA /
        (1 + theta * H0Yij_HA) - sum(as.numeric((Yij_cal_HA < cens_time)) * shape.W / scale_W),
      -(1 / theta + Di_) * theta * dH0_epsiYij_HA /
        (1 + theta * H0Yij_HA) - sum(as.numeric((Yij_HA < cens_time)) * shape.W / scale_W)
    )

    nu_fact_HA <- ifelse((tolower(data.type) == "rec_event"),
      sum(c(as.numeric((Yij_cal_HA < cens_time)) * (dh0_nu(t = Yij_HA, y = scale_W, p = shape.W) / h0(t = Yij_HA, y = scale_W, p = shape.W))), na.rm = TRUE),
      sum(as.numeric((Yij_HA < cens_time)) * (dh0_nu(t = Yij_HA, y = scale_W, p = shape.W) / h0(t = Yij_HA, y = scale_W, p = shape.W)))
    )
    score_nu_HA <- -(1 / theta + Di) * theta * dH0_nuYij_HA / (1 + theta * H0Yij_HA) + nu_fact_HA

    if (npar == 4) {
      score_vector_H0 <- c(score_beta_H0, score_theta_H0, score_epsi_H0, score_nu_H0)
      score_vector_HA <- c(score_beta_HA, score_theta_HA, score_epsi_HA, score_nu_HA)
    }

    U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
    U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)


    sum_score2_H0 <- sum_score2_H0 + U_score2_H0
    sum_score2_HA <- sum_score2_HA + U_score2_HA

    Di_H0 <- Di_H0 + Di
    Di_HA <- Di_HA + Di_
    sum_ni <- sum_ni + ni_

    if (Di == 0) {
      cens_E <- cens_E + 1
    }
    if (Di_ == 0) {
      cens_E_ <- cens_E_ + 1
    }
  }
  E_score2_H0 <- sum_score2_H0 / samples.mc
  E_score2_HA <- sum_score2_HA / samples.mc
  E_Di_H0 <- Di_H0 / sum_ni
  E_Di_HA <- Di_HA / sum_ni

  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }

  C_trt <- matrix(data = c(1, 0, 0, 0), nrow = 1, ncol = npar)

  if (statistic == "Wald") {
    ncp0 <- (Gc) * t(C_trt %*% (c(beta.H0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(beta.H0, 0, 0, 0)))
    ncpA <- (Gc) * t(C_trt %*% (c(beta.HA, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0)))

    rr <- nrow(C_trt)
  }

  res <- list(
    estimated.power     = 1 - pchisq(qchisq(1 - alpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA),
    events    = c(Gc * Di_H0 / samples.mc, Gc * Di_HA / samples.mc)
    # censoring = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}

Nsn_SFM <- function(
    power = 0.8, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1,
    beta.H0 = 0, beta.HA = log(0.75),
    shape.W = 1, theta = 0.25, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo",
    statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(power)) {
    stop("Power is missing")
  }
  if ((power <= 0) | (power > 0.99)) {
    stop("power must be between ]0,1[")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'unif' or 'pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is fixed ")
  }
  if (data.type == "grouped") {
    if ((length(ni) == 1) & (ni[1] <= 0)) {
      stop("ni must be strictly positive")
    }
    if ((length(ni) == 2) & ((ni[1] <= 0) | (ni[2] <= 0))) {
      stop("ni must be strictly positive")
    }
  } else {
    if ((length(ni) == 1) & (ni[1] < 0)) {
      stop("ni must be  positive")
    }
    if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
      stop("ni must be  positive")
    }
  }
  if ((!tolower(data.type) %in% c("rec_event", "grouped"))) {
    stop("data.type must be either 'grouped' or 'rec_event' ")
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be either 'gap' or 'calendar'")
  }
  if (tolower(data.type) == "grouped" & tolower(timescale) != "gap") {
    stop("For grouped data structure, one must use a gap timescale")
  }
  if (tolower(data.type) == "grouped" & tolower(ni.type) != "max") {
    stop("For grouped data structure, ni.type must be max (to ensure a fixed number of subjects per group)")
  }
  if ((!tolower(cens.type) %in% c ("unif", "uniform", "u", "exp", "expo", "exponential"))) {
    stop("cens.type must be either 'Uniform' or 'Exponential ")
  }
  if ((tolower(cens.type) %in% c("u", "unif", "uniform")) & length(cens.par) != 2) {
    stop("cens.par must be a vector of length two if cens.type is uniform ")
  }
  if ((tolower(cens.type) %in% c("expo", "exp", "exponential")) & length(cens.par) != 1) {
    stop("cens.par must be a scalar if cens.type is exponential ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }

  P_E <- 1 / (1 / ratio + 1)

  theta <- theta

  beta.HA <- beta.HA
  beta.H0 <- beta.H0

  npar <- 4
  if (test.type == "1-sided") {
    alpha <- typeIerror
  }
  if (test.type == "2-sided") {
    alpha <- typeIerror / 2
  }

  scale_W <- median.H0 / log(2)^(1 / shape.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }

  cumsumrev <- function(x) {
    output <- rep(NA, length(x))
    ll <- length(x)
    if (ll > 1) {
      x1 <- x[1:(ll - 1)]
      x2 <- x[2:ll]
      output <- c(x[1], x2 - x1)
    }
    if (ll == 1) {
      output <- x
    }
    return(output)
  }

  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_ni <- 0
  Di_H0 <- Di_HA <- 0

  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "fixed"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }
    if (ni_ > 0) {
      if (data.type == "rec_event") {
        Xij <- rep(sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E)), ni_)
        Accrual <- rep(runif(n = 1, min = 0, max = Acc.Dur), ni_)
      } else {
        Xij <- sample(x = c(0, 1), size = ni_, prob = c(1 - P_E, P_E), replace = TRUE)
        Accrual <- runif(n = ni_, min = 0, max = Acc.Dur)
      }

      U <- runif(n = ni_, 0, 1)
      Vi <- rgamma(1, shape = theta^(-1), scale = theta)

      Tij_H0 <- H0_1(s = -log(U) * exp(-(beta.H0 * Xij)) / Vi, y = scale_W, p = shape.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(beta.HA * Xij)) / Vi, y = scale_W, p = shape.W)

      if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
        if (data.type == "rec_event") {
          deathtime <- runif(1, cens.par[1], cens.par[2])
        } else {
          deathtime <- runif(ni_, cens.par[1], cens.par[2])
        }
      } else {
        death_rate <- log(2) / cens.par
        if (data.type == "rec_event") {
          deathtime <- rexp(1, death_rate)
        } else {
          deathtime <- rexp(ni_, death_rate)
        }
      }

      if (tolower(FUP.type) %in% c("f", "fixed")) {
        cens_time <- pmin(FUP, deathtime)
      }
      if (tolower(FUP.type) %in%  c("uptoend", "end")) {
        cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
      }
    } else {
      Xij <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
      Accrual <- runif(n = 1, min = 0, max = Acc.Dur)

      if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
        deathtime <- runif(1, cens.par[1], cens.par[2])
      } else {
        death_rate <- log(2) / cens.par
        deathtime <- rexp(1, death_rate)
      }
      if (FUP.type %in% c("Fixed", "fixed", "f")) {
        cens_time <- pmin(FUP, deathtime)
      }
      if (tolower(FUP.type) %in%  c("uptoend", "end")) {
        cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
      }
      Tij_H0 <- cens_time
      Tij_HA <- cens_time
    }


    if (data.type == "rec_event") {
      Yij_cal_H0 <- Yij_H0 <- pmin(cumsum(Tij_H0), cens_time)
      Yij_cal_HA <- Yij_HA <- pmin(cumsum(Tij_HA), cens_time)

      if (tolower(timescale) == "gap") {
        Yij_H0 <- cumsumrev(Yij_cal_H0)
        Yij_HA <- cumsumrev(Yij_cal_HA)
      } else {
        Yij_H0 <- Yij_cal_H0
        Yij_HA <- Yij_cal_HA
      }
      Di <- sum(as.numeric((Yij_cal_H0 < cens_time)))
      Di_ <- sum(as.numeric((Yij_cal_HA < cens_time)))
      DiXij <- sum(as.numeric((Yij_cal_H0 < cens_time)) * Xij)
      DiXij_ <- sum(as.numeric((Yij_cal_HA < cens_time)) * Xij)
    } else {
      Yij_H0 <- pmin(Tij_H0, cens_time)
      Yij_HA <- pmin(Tij_HA, cens_time)
      Di <- sum(as.numeric(Yij_H0 < cens_time))
      Di_ <- sum(as.numeric(Yij_HA < cens_time))
      DiXij <- sum(as.numeric(Yij_H0 < cens_time) * Xij)
      DiXij_ <- sum(as.numeric(Yij_HA < cens_time) * Xij)
    }

    H0Yij_H0 <- sum(c(H0(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))))
    dH0_epsiYij_H0 <- sum(c(dH0_epsi(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))))
    dH0_nuYij_H0 <- sum(c(dH0_nu(t = Yij_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xij))), na.rm = TRUE)

    H0Yij_HA <- sum(c(H0(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))))
    dH0_epsiYij_HA <- sum(c(dH0_epsi(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))))
    dH0_nuYij_HA <- sum(c(dH0_nu(t = Yij_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xij))), na.rm = TRUE)


    ln <- c(1:(Di))
    ln_ <- c(1:(Di_))

    score_beta_H0 <- DiXij - (1 / theta + Di) * (theta * sum(H0(t = Yij_H0, y = scale_W, p = shape.W) * Xij *
      exp(beta.H0 * Xij))) / (1 + theta * (H0Yij_H0))


    score_theta_H0 <- log(1 + theta * H0Yij_H0) / theta^2 - (1 / theta + Di) * H0Yij_H0 / (1 + theta * H0Yij_H0) +
      ifelse(Di > 0, sum((Di - ln) / (1 + theta * (Di - ln))), 0)

    score_epsi_H0 <- ifelse((tolower(data.type) == "rec_event"),
      -(1 / theta + Di) * theta * dH0_epsiYij_H0 /
        (1 + theta * H0Yij_H0) - sum(as.numeric((Yij_cal_H0 < cens_time)) * shape.W / scale_W),
      -(1 / theta + Di) * theta * dH0_epsiYij_H0 /
        (1 + theta * H0Yij_H0) - sum(as.numeric((Yij_H0 < cens_time)) * shape.W / scale_W)
    )

    nu_fact_H0 <- ifelse((tolower(data.type) == "rec_event"),
      sum(c(as.numeric((Yij_cal_H0 < cens_time)) * (dh0_nu(t = Yij_H0, y = scale_W, p = shape.W) / h0(t = Yij_H0, y = scale_W, p = shape.W))), na.rm = TRUE),
      sum(as.numeric((Yij_H0 < cens_time)) * (dh0_nu(t = Yij_H0, y = scale_W, p = shape.W) / h0(t = Yij_H0, y = scale_W, p = shape.W)))
    )
    score_nu_H0 <- -(1 / theta + Di) * theta * dH0_nuYij_H0 / (1 + theta * H0Yij_H0) + nu_fact_H0

    score_beta_HA <- DiXij_ - (1 / theta + Di_) * (theta * sum(H0(t = Yij_HA, y = scale_W, p = shape.W) * Xij *
      exp(beta.HA * Xij))) / (1 + theta * (H0Yij_HA))

    score_theta_HA <- log(1 + theta * H0Yij_HA) / theta^2 - (1 / theta + Di_) * H0Yij_HA / (1 + theta * H0Yij_HA) +
      ifelse(Di_ > 0, sum((Di_ - ln_) / (1 + theta * (Di_ - ln_))), 0)

    score_epsi_HA <- ifelse((tolower(data.type) == "rec_event"),
      -(1 / theta + Di_) * theta * dH0_epsiYij_HA /
        (1 + theta * H0Yij_HA) - sum(as.numeric((Yij_cal_HA < cens_time)) * shape.W / scale_W),
      -(1 / theta + Di_) * theta * dH0_epsiYij_HA /
        (1 + theta * H0Yij_HA) - sum(as.numeric((Yij_HA < cens_time)) * shape.W / scale_W)
    )

    nu_fact_HA <- ifelse((tolower(data.type) == "rec_event"),
      sum(c(as.numeric((Yij_cal_HA < cens_time)) * (dh0_nu(t = Yij_HA, y = scale_W, p = shape.W) / h0(t = Yij_HA, y = scale_W, p = shape.W))), na.rm = TRUE),
      sum(as.numeric((Yij_HA < cens_time)) * (dh0_nu(t = Yij_HA, y = scale_W, p = shape.W) / h0(t = Yij_HA, y = scale_W, p = shape.W)))
    )
    score_nu_HA <- -(1 / theta + Di) * theta * dH0_nuYij_HA / (1 + theta * H0Yij_HA) + nu_fact_HA

    if (npar == 4) {
      score_vector_H0 <- c(score_beta_H0, score_theta_H0, score_epsi_H0, score_nu_H0)
      score_vector_HA <- c(score_beta_HA, score_theta_HA, score_epsi_HA, score_nu_HA)
    }

    U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
    U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)


    sum_score2_H0 <- sum_score2_H0 + U_score2_H0
    sum_score2_HA <- sum_score2_HA + U_score2_HA

    if (ni_ > 0) {
      Di_H0 <- Di_H0 + Di
      Di_HA <- Di_HA + Di_
      sum_ni <- sum_ni + ni_
    }
  }
  E_score2_H0 <- sum_score2_H0 / samples.mc
  E_score2_HA <- sum_score2_HA / samples.mc
  E_Di_H0 <- Di_H0 / sum_ni
  E_Di_HA <- Di_HA / sum_ni
  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }

  C_trt <- matrix(data = c(1, 0, 0, 0), nrow = 1, ncol = npar)
  rr <- nrow(C_trt)

  ncp <- .chi2ncp(alpha * 2, power, df = rr)
  CIC <- t(C_trt %*% (c(beta.HA, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0)))
  if (statistic == "Wald") {
    Npts <- ceiling(ncp * solve(CIC))
    ncp0 <- (Npts) * t(C_trt %*% (c(beta.H0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(beta.H0, 0, 0, 0)))
    ncpA <- (Npts) * t(C_trt %*% (c(beta.HA, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0)))
    estimated.power <- 1 - pchisq(qchisq(1 - alpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA)
  }

  res <- list(
    Groups = Npts,
    estimated.power = estimated.power,
    events = c(Npts * Di_H0 / samples.mc, Npts * Di_HA / samples.mc)
    # censoring  = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}





Power_NFM <- function(
    Groups = 80, ni = 8, ni.type = "max", kij = 15, kij.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1, beta.H0 = 0, beta.HA = log(0.75),
    shape.W = 1, theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo", statistic = "Wald",
    typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(Groups)) {
    stop("Groups is missing")
  }
  if (Groups <= 0) {
    stop("Groups must be strictly positive")
  }
  if (tolower(data.type) == "grouped" & (tolower(ni.type) != "max" | tolower(kij.type) != "max")) {
    stop("For grouped data structure, both ni.type and kij.type must be max (to ensure a
    fixed number of subgroups per group, and a fixed number of subjects per subgroups)")
  }
  if (tolower(data.type) == "rec_event1" & (tolower(ni.type) != "max")) {
    stop("Under 'rec_event1', ni corresponds to a number of subjects. In order to ensure a
    fixed number of subjects per groups, only 'max' option is valid.")
  }
  if ((data.type == "rec_event2" & tolower(ni.type) != "max")) {
    stop("Under 'rec_event2', ni corresponds to a number of distinct recurrent events type. 
    Hence, only 'max' is valid.")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'unif'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is fixed ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & ((ni[1] <= 0) | (ni[2] <= 0))) {
    stop("ni must be strictly positive")
  }
  if ((!tolower(kij.type) %in% c("maximum", "max", "m", "unif", "uniform", "u", "p", "pois", "poisson"))) {
    stop("kij.type must be either 'max' or 'unif' or 'pois'")
  }
  if ((tolower(kij.type) %in% c("u", "unif", "uniform")) & length(kij) != 2) {
    stop("kij must be a vector of length two if kij.type is uniform")
  }
  if ((tolower(kij.type) %in% c("fixed", "Fixed", "F", "f")) & length(kij) != 1) {
    stop("kij must be a scalar if kij.type is fixed")
  }
  if ((length(kij) == 1) & (kij[1] <= 0)) {
    stop("kij must be strictly positive")
  }
  if ((length(kij) == 2) & ((kij[1] <= 0) | (kij[2] <= 0))) {
    stop("kij must be strictly positive")
  }
  if ((!tolower(data.type) %in% c("grouped", "rec_event1", "rec_event2"))) {
    stop("data.type must be either 'grouped', 'rec_event1' or 'rec_event2' ")
  }
  if ((data.type %in% c("grouped")) & (min(kij) <= 0)) {
    stop("kij must be strictly positive when data.type is grouped")
  }
  if (timescale %in% c("Calendar", "calendar") & data.type != "rec_event1") {
    stop("Calendar timescale is only possible when one type of recurrent event is considered")
  }
  if ((!tolower(cens.type) %in% c ("unif", "uniform", "u", "exp", "expo", "exponential"))) {
    stop("cens.type must be either 'Uniform' or 'Exponential")
  }
  if ((tolower(cens.type) %in% c("u", "unif", "uniform")) & length(cens.par) != 2) {
    stop("cens.par must be a vector of length two if cens.type is uniform ")
  }
  if ((tolower(cens.type) %in% c("expo", "exp", "exponential")) & length(cens.par) != 1) {
    stop("cens.par must be a scalar if cens.type is exponential ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }


  Gc <- Groups
  P_E <- 1 / (1 / ratio + 1)

  eta <- eta
  theta <- theta
  beta.HA <- beta.HA
  beta.H0 <- beta.H0

  npar <- 5
  if (test.type == "1-sided") {
    alpha <- typeIerror
  }
  if (test.type == "2-sided") {
    alpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scale_W <- median.H0 / log(2)^(1 / shape.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }
  cumsumrev <- function(x, id) {
    ind <- unique(id)
    output <- rep(NA, length(x))
    for (ii in c(1:length(ind))) {
      indice <- which(ii == id)
      xind <- x[indice]
      ll <- length(xind)
      if (ll > 1) {
        x1 <- xind[1:(ll - 1)]
        x2 <- xind[2:ll]
        output[indice] <- c(xind[1], x2 - x1)
      }
      if (ll == 1) {
        output[indice] <- xind
      }
    }
    return(output)
  }
  cumsum_m <- function(x, id) {
    ind <- unique(id)
    output <- rep(NA, length(x))
    for (ii in c(1:length(ind))) {
      indice <- which(ii == id)
      xind <- x[indice]
      output[indice] <- cumsum(xind)
    }
    return(output)
  }
  G_THETA <- function(vi, Di, eta, theta, H00, Dij) {
    return(vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta)) / prod(((1 + theta * vi * H00)^{
        1 / theta + Dij
      })))
  }
  dbeta_G_THETA <- function(vi, eta, theta, Di, Dij, H00, H00_X) {
    d1 <- (1 / theta + Dij) * (1 + theta * vi * H00)^{
      1 / theta + Dij - 1
    } * theta * vi * H00_X
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }
    dbeta_1 <- prod(f1) * sum(d1 / f1)

    return(-vi^{
      Di + 1 / eta - 1
    } * exp(vi * (1 - 1 / eta)) * dbeta_1 / prod(f1)^2)
  }

  deta_G_THETA <- function(vi, eta, theta, Di, Dij, H00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta)) * (vi - log(vi)) / eta^2
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }
    return(d1 / prod(f1))
  }
  dtheta_G_THETA <- function(vi, eta, theta, Di, Dij, H00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta))
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }

    p1 <- f1 * ((-1 / theta^2) * log(1 + theta * vi * H00) + {
      1 / theta + Dij
    } * vi * H00 / (1 + theta * vi * H00))
    dtheta_1 <- prod(f1) * sum(p1 / f1)

    return(-d1 * dtheta_1 / prod(f1)^2)
  }


  dh0_G_THETA <- function(vi, eta, theta, Di, Dij, H00, dH00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta))
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }

    p1 <- (Dij + 1 / theta) * (1 + theta * vi * H00)^{
      Dij + 1 / theta - 1
    } * vi * theta * dH00
    dh0_1 <- prod(f1) * sum(p1 / f1)

    return(-d1 * dh0_1 / prod(f1)^2)
  }


  sum_score2_H0 <- sum_score2_HA <- U_score2_H0 <- U_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- sumNA <- 0
  Di_H0 <- Di_HA <- 0


  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "fixed"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }


    if (tolower(kij.type) %in% c("u", "unif", "uniform")) {
      if (kij[1] == kij[2]) {
        kij.type <- "fixed"
        tmp <- kij[1]
        kij <- tmp
        rm(tmp)
      } else {
        kij_ <- round(sample(kij[1]:kij[2], ni_, replace = TRUE))
      }
    }
    if (tolower(kij.type) %in% c("Maximum", "maximum", "Max", "max")) {
      kij_ <- rep(round(kij), ni_)
    }

    if (tolower(kij.type) %in% c("p", "pois", "poisson")) {
      kij_ <- rpois(n = ni_, lambda = kij)
    }

    nb.zerok <- sum(kij_ == 0)
    if (data.type == "rec_event2") {
      if (nb.zerok > 0) {
        where.zerok <- which(kij_ == 0)
        kij_[where.zerok] <- 1
      }
      Xijk <- rep(sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E), replace = TRUE), sum(kij_))
      Accrual <- rep(runif(n = 1, min = 0, max = Acc.Dur), sum(kij_))
    }
    if (data.type == "rec_event1") {
      if (nb.zerok > 0) {
        where.zerok <- which(kij_ == 0)
        kij_[where.zerok] <- 1
      }
      Xijk <- rep(sample(x = c(0, 1), size = ni_, prob = c(1 - P_E, P_E), replace = TRUE), kij_)
      Accrual <- rep(runif(n = ni_, min = 0, max = Acc.Dur), kij_)
    }
    if (data.type == "grouped") {
      if (nb.zerok > 0) {
        ni_ <- ni_ - nb.zerok
        kij_ <- kij_[kij_ != 0]
      }
      if (ni_ <= 0) {
        stop("Sample size too small, empty subgroups ...")
      }
      Xijk <- sample(x = c(0, 1), size = sum(kij_), prob = c(1 - P_E, P_E), replace = TRUE)
      Accrual <- runif(n = sum(kij_), min = 0, max = Acc.Dur)
    }

    ID_subgroup <- rep(c(1:ni_), kij_)

    U <- runif(n = sum(kij_), 0, 1) ##

    vi <- rgamma(1, shape = eta^(-1), scale = eta)
    wij <- rep(rgamma(ni_, shape = theta^(-1), scale = theta), kij_)

    if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
      if (data.type == "rec_event2") {
        deathtime <- runif(1, cens.par[1], cens.par[2])
      }
      if (data.type == "rec_event1") {
        deathtime <- rep(runif(ni_, cens.par[1], cens.par[2]), kij_)
      }
      if (data.type == "grouped") {
        deathtime <- runif(sum(kij_), cens.par[1], cens.par[2])
      }
    } else {
      death_rate <- log(2) / cens.par
      if (data.type == "rec_event2") {
        deathtime <- rexp(1, death_rate)
      }
      if (data.type == "rec_event1") {
        deathtime <- rep(rexp(ni_, death_rate), kij_)
      }
      if (data.type == "grouped") {
        deathtime <- runif(sum(kij_), death_rate)
      }
    }

    if (tolower(FUP.type) %in% c("f", "fixed")) {
      cens_time <- pmin(FUP, deathtime)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
    }

    Tijk_H0 <- H0_1(s = -log(U) * exp(-(beta.H0 * Xijk)) / (vi * wij), y = scale_W, p = shape.W)
    Tijk_HA <- H0_1(s = -log(U) * exp(-(beta.HA * Xijk)) / (vi * wij), y = scale_W, p = shape.W)

    ind <- unique(ID_subgroup)
    ind_cum_sub <- sapply(ind, function(i) which(ID_subgroup == i)[length(which(ID_subgroup == i))])

    if ((data.type == "rec_event2") || (data.type == "rec_event1")) {
      if (nb.zerok > 0) {
        Tijk_H0[ind_cum_sub[where.zerok]] <- cens_time[ind_cum_sub[where.zerok]]
        Tijk_HA[ind_cum_sub[where.zerok]] <- cens_time[ind_cum_sub[where.zerok]]
      }
      Yijk_cal_H0 <- pmin(cumsum_m(Tijk_H0, ID_subgroup), cens_time)
      Yijk_cal_HA <- pmin(cumsum_m(Tijk_HA, ID_subgroup), cens_time)

      if (timescale %in% c("calendar", "Calendar") & data.type == "rec_event1") {
        Yijk_H0 <- Yijk_cal_H0
        Yijk_HA <- Yijk_cal_HA
      } else {
        Yijk_H0 <- cumsumrev(Yijk_cal_H0, ID_subgroup)
        Yijk_HA <- cumsumrev(Yijk_cal_HA, ID_subgroup)
      }
      deltaijk <- as.numeric((Yijk_cal_H0 < cens_time))
      deltaijk_ <- as.numeric((Yijk_cal_HA < cens_time))

      Di <- sum(deltaijk)
      Di_ <- sum(deltaijk_)
      Dij <- ave(deltaijk, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      Dij_ <- ave(deltaijk_, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      DiXijk <- sum(deltaijk * Xijk)
      DiXijk_ <- sum(deltaijk_ * Xijk)
    } else {
      Yijk_H0 <- pmin(Tijk_H0, cens_time)
      Yijk_HA <- pmin(Tijk_HA, cens_time)

      deltaijk <- as.numeric((Yijk_H0 < cens_time))
      deltaijk_ <- as.numeric((Yijk_HA < cens_time))

      Di <- sum(deltaijk)
      Di_ <- sum(deltaijk_)
      Dij <- ave(deltaijk, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      Dij_ <- ave(deltaijk_, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      DiXijk <- sum(deltaijk * Xijk)
      DiXijk_ <- sum(deltaijk_ * Xijk)
    }


    H0Yijk_H0 <- ave(H0(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    H0Yijk_HA <- ave(H0(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    dH0_epsiYijk_H0 <- ave(dH0_epsi(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    dH0_epsiYijk_HA <- ave(dH0_epsi(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    dH0_nuYijk_H0 <- ave(dH0_nu(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = function(x) sum(x, na.rm = T))[ind_cum_sub]
    dH0_nuYijk_HA <- ave(dH0_nu(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = function(x) sum(x, na.rm = T))[ind_cum_sub]

    H0Yijk_H0_X <- ave(H0(t = Yijk_H0, y = scale_W, p = shape.W) * Xijk * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    H0Yijk_HA_X <- ave(H0(t = Yijk_HA, y = scale_W, p = shape.W) * Xijk * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))

    dbeta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbeta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, H00_X = H0Yijk_H0_X, Dij = Dij)
    }))
    dbeta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbeta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, H00_X = H0Yijk_HA_X, Dij = Dij_)
    }))
    score_beta_H0 <- DiXijk + dbeta_G_THETA_val_H0 / G_THETA_val_H0
    score_beta_HA <- DiXijk_ + dbeta_G_THETA_val_HA / G_THETA_val_HA

    deta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    deta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))

    score_eta_H0 <- digamma(1 / eta) / eta^2 + (log(eta) - 1) / eta^2 + deta_G_THETA_val_H0 / G_THETA_val_H0
    score_eta_HA <- digamma(1 / eta) / eta^2 + (log(eta) - 1) / eta^2 + deta_G_THETA_val_HA / G_THETA_val_HA

    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))

    indic <- sum(sapply(c(1:length(kij_)), function(jj) {
      if (Dij[jj] > 0) {
        ln <- c(1:(sum(Dij[jj])))
        return(sum((Dij[jj] - ln) / (1 + theta * (Dij[jj] - ln))))
      } else {
        0
      }
    }))
    indic_ <- sum(sapply(c(1:length(kij_)), function(jj) {
      if (Dij_[jj] > 0) {
        ln_ <- c(1:(sum(Dij_[jj])))
        return(sum((Dij_[jj] - ln_) / (1 + theta * (Dij_[jj] - ln_))))
      } else {
        0
      }
    }))
    score_theta_H0 <- indic + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- indic_ + dtheta_G_THETA_val_HA / G_THETA_val_HA

    depsi_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, dH00 = dH0_epsiYijk_H0, Dij = Dij)
    }))
    depsi_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, dH00 = dH0_epsiYijk_HA, Dij = Dij_)
    }))


    score_epsi_H0 <- -Di * shape.W / scale_W + depsi_G_THETA_val_H0 / G_THETA_val_H0
    score_epsi_HA <- -Di_ * shape.W / scale_W + depsi_G_THETA_val_HA / G_THETA_val_HA

    dnu_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, dH00 = dH0_nuYijk_H0, Dij = Dij)
    }))
    dnu_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, dH00 = dH0_nuYijk_HA, Dij = Dij_)
    }))

    score_nu_H0 <- sum(deltaijk * dh0_nu(Yijk_H0, scale_W, shape.W) / h0(Yijk_H0, scale_W, shape.W), na.rm = TRUE) + dnu_G_THETA_val_H0 / G_THETA_val_H0
    score_nu_HA <- sum(deltaijk_ * dh0_nu(Yijk_HA, scale_W, shape.W) / h0(Yijk_HA, scale_W, shape.W), na.rm = TRUE) + dnu_G_THETA_val_HA / G_THETA_val_HA

    if (npar == 5) {
      score_vector_H0 <- c(score_beta_H0, score_eta_H0, score_theta_H0, score_epsi_H0, score_nu_H0)
      score_vector_HA <- c(score_beta_HA, score_eta_HA, score_theta_HA, score_epsi_HA, score_nu_HA)
    }
    if (sum(is.na(c(score_vector_H0, score_vector_HA))) == 0) {
      U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
      U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)
      sum_score2_H0 <- sum_score2_H0 + U_score2_H0
      sum_score2_HA <- sum_score2_HA + U_score2_HA

      sum_Di_H0 <- sum_Di_H0 + Di / sum(kij_)
      Di_H0 <- Di_H0 + Di
      sum_Di_HA <- sum_Di_HA + Di_ / sum(kij_)
      Di_HA <- Di_HA + Di_
    } else {
      sumNA <- sumNA + 1
    }
  }
  E_score2_H0 <- sum_score2_H0 / (samples.mc - sumNA)
  E_score2_HA <- sum_score2_HA / (samples.mc - sumNA)
  E_Di_H0 <- sum_Di_H0 / (samples.mc - sumNA)
  E_Di_HA <- sum_Di_HA / (samples.mc - sumNA)


  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }


  C_trt <- matrix(data = c(1, 0, 0, 0, 0), nrow = 1, ncol = npar)

  cat("\n")

  if (statistic == "Wald") {
    ncp0 <- (Gc) * t(C_trt %*% (c(beta.H0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(beta.H0, 0, 0, 0, 0)))
    ncpA <- (Gc) * t(C_trt %*% (c(beta.HA, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0, 0)))
    rr <- nrow(C_trt)
  }

  res <- list(
    estimated.power = 1 - pchisq(qchisq(1 - alpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA),
    events          = c(Gc * Di_H0 / (samples.mc - sumNA), Gc * Di_HA / (samples.mc - sumNA))
    # censoring     = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}


Nsn_NFM <- function(
    power = 0.8, ni = 8, ni.type = "max", kij = 15, kij.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", median.H0 = 1, beta.H0 = 0, beta.HA = log(0.75),
    shape.W = 1, theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", data.type = "grouped", cens.par = 5, cens.type = "Expo", statistic = "Wald",
    typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(power)) {
    stop("Power is missing")
  }
  if (tolower(data.type) == "grouped" & (tolower(ni.type) != "max" | tolower(kij.type) != "max")) {
    stop("For grouped data structure, both ni.type and kij.type must be max (to ensure a
    fixed number of subgroups per group, and a fixed number of subjects per subgroups)")
  }
  if (tolower(data.type) == "rec_event1" & (tolower(ni.type) != "max")) {
    stop("Under 'rec_event1', ni corresponds to a number of subjects. In order to ensure a
    fixed number of subjects per groups, only 'max' option is valid.")
  }
  if ((tolower(data.type) == "rec_event2" & tolower(ni.type) != "max")) {
    stop("Under 'rec_event2', ni corresponds to a number of distinct recurrent events type. 
    Hence, only 'max' is valid.")
  }
  if ((power <= 0) | (power > 0.99)) {
    stop("power must be between ]0,1[")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'Unif' or 'Pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is fixed ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & ((ni[1] <= 0) | (ni[2] <= 0))) {
    stop("ni must be strictly positive")
  }
  if ((!tolower(kij.type) %in% c("maximum", "max", "m", "unif", "uniform", "u", "poisson", "pois", "p"))) {
    stop("kij.type must be either 'max' or 'Unif' or 'Pois'")
  }
  if ((tolower(kij.type) %in% c("u", "unif", "uniform")) & length(kij) != 2) {
    stop("kij must be a vector of length two if kij.type is uniform ")
  }
  if ((tolower(kij.type) %in% c("fixed", "f")) & length(kij) != 1) {
    stop("kij must be a scalar if kij.type is fixed ")
  }
  if ((length(kij) == 1) & (kij[1] <= 0)) {
    stop("kij must be strictly positive")
  }
  if ((length(kij) == 2) & ((kij[1] <= 0) | (kij[2] <= 0))) {
    stop("kij must be strictly positive")
  }
  if ((!tolower(data.type) %in% c("grouped", "rec_event1", "rec_event2"))) {
    stop("data.type must be either 'grouped', 'rec_event1' or 'rec_event2' ")
  }
  if ((data.type %in% c("grouped")) & (min(kij) <= 0)) {
    stop("kij must be strictly positive when data.type is grouped")
  }
  if (timescale %in% c("Calendar", "calendar") & data.type == "rec_event1") {
    stop("Calendar timescale is only possible when one type of recurrent event is considered")
  }
  if ((!tolower(cens.type) %in% c ("unif", "uniform", "u", "exp", "expo", "exponential"))) {
    stop("cens.type must be either 'Uniform' or 'Exponential ")
  }
  if ((tolower(cens.type) %in% c("u", "unif", "uniform")) & length(cens.par) != 2) {
    stop("cens.par must be a vector of length two if cens.type is uniform ")
  }
  if ((tolower(cens.type) %in% c("expo", "exp", "exponential")) & length(cens.par) != 1) {
    stop("cens.par must be a scalar if cens.type is exponential ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }

  P_E <- 1 / (1 / ratio + 1)
  eta <- eta
  theta <- theta
  beta.HA <- beta.HA
  beta.H0 <- beta.H0

  npar <- 5
  if (test.type == "1-sided") {
    alpha <- typeIerror
  }
  if (test.type == "2-sided") {
    alpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scale_W <- median.H0 / log(2)^(1 / shape.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }
  cumsumrev <- function(x, id) {
    ind <- unique(id)
    output <- rep(NA, length(x))
    for (ii in c(1:length(ind))) {
      indice <- which(ii == id)
      xind <- x[indice]
      ll <- length(xind)
      if (ll > 1) {
        x1 <- xind[1:(ll - 1)]
        x2 <- xind[2:ll]
        output[indice] <- c(xind[1], x2 - x1)
      }
      if (ll == 1) {
        output[indice] <- xind
      }
    }
    return(output)
  }
  cumsum_m <- function(x, id) {
    ind <- unique(id)
    output <- rep(NA, length(x))
    for (ii in c(1:length(ind))) {
      indice <- which(ii == id)
      xind <- x[indice]
      output[indice] <- cumsum(xind)
    }
    return(output)
  }
  G_THETA <- function(vi, Di, eta, theta, H00, Dij) {
    return(vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta)) / prod(((1 + theta * vi * H00)^{
        1 / theta + Dij
      })))
  }
  dbeta_G_THETA <- function(vi, eta, theta, Di, Dij, H00, H00_X) {
    d1 <- (1 / theta + Dij) * (1 + theta * vi * H00)^{
      1 / theta + Dij - 1
    } * theta * vi * H00_X
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }

    dbeta_1 <- prod(f1) * sum(d1 / f1)

    return(-vi^{
      Di + 1 / eta - 1
    } * exp(vi * (1 - 1 / eta)) * dbeta_1 / prod(f1)^2)
  }

  deta_G_THETA <- function(vi, eta, theta, Di, Dij, H00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta)) * (vi - log(vi)) / eta^2
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }
    return(d1 / prod(f1))
  }
  dtheta_G_THETA <- function(vi, eta, theta, Di, Dij, H00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta))
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }


    p1 <- f1 * ((-1 / theta^2) * log(1 + theta * vi * H00) + {
      1 / theta + Dij
    } * vi * H00 / (1 + theta * vi * H00))

    dtheta_1 <- prod(f1) * sum(p1 / f1)

    return(-d1 * dtheta_1 / prod(f1)^2)
  }


  dh0_G_THETA <- function(vi, eta, theta, Di, Dij, H00, dH00) {
    d1 <- vi^
      {
        Di + 1 / eta - 1
      } * exp(vi * (1 - 1 / eta))
    f1 <- (1 + theta * vi * H00)^{
      1 / theta + Dij
    }


    p1 <- (Dij + 1 / theta) * (1 + theta * vi * H00)^{
      Dij + 1 / theta - 1
    } * vi * theta * dH00

    dh0_1 <- prod(f1) * sum(p1 / f1)

    return(-d1 * dh0_1 / prod(f1)^2)
  }


  sum_score2_H0 <- sum_score2_HA <- U_score2_H0 <- U_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- sumNA <- 0
  Di_H0 <- Di_HA <- 0


  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "fixed"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }


    if (tolower(kij.type) %in% c("u", "unif", "uniform")) {
      if (kij[1] == kij[2]) {
        kij.type <- "fixed"
        tmp <- kij[1]
        kij <- tmp
        rm(tmp)
      } else {
        kij_ <- round(sample(kij[1]:kij[2], ni_, replace = TRUE))
      }
    }
    if (tolower(kij.type) %in% c("Maximum", "maximum", "Max", "max")) {
      kij_ <- rep(round(kij), ni_)
    }

    if (tolower(kij.type) %in% c("p", "pois", "poisson")) {
      kij_ <- rpois(n = ni_, lambda = kij)
    }

    nb.zerok <- sum(kij_ == 0)
    if (data.type == "rec_event2") {
      if (nb.zerok > 0) {
        where.zerok <- which(kij_ == 0)
        kij_[where.zerok] <- 1
      }
      Xijk <- rep(sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E), replace = TRUE), sum(kij_))
      Accrual <- rep(runif(n = 1, min = 0, max = Acc.Dur), sum(kij_))
    }
    if (data.type == "rec_event1") {
      if (nb.zerok > 0) {
        where.zerok <- which(kij_ == 0)
        kij_[where.zerok] <- 1
      }
      Xijk <- rep(sample(x = c(0, 1), size = ni_, prob = c(1 - P_E, P_E), replace = TRUE), kij_)
      Accrual <- rep(runif(n = ni_, min = 0, max = Acc.Dur), kij_)
    }
    if (data.type == "grouped") {
      if (nb.zerok > 0) {
        ni_ <- ni_ - nb.zerok
        kij_ <- kij_[kij_ != 0]
      }
      if (ni_ <= 0) {
        stop("Sample size too small, empty subgroups ...")
      }
      Xijk <- sample(x = c(0, 1), size = sum(kij_), prob = c(1 - P_E, P_E), replace = TRUE)
      Accrual <- runif(n = sum(kij_), min = 0, max = Acc.Dur)
    }

    ID_subgroup <- rep(c(1:ni_), kij_)

    U <- runif(n = sum(kij_), 0, 1) ##

    vi <- rgamma(1, shape = eta^(-1), scale = eta)
    wij <- rep(rgamma(ni_, shape = theta^(-1), scale = theta), kij_)

    if (tolower(cens.type) %in% c("u", "unif", "uniform")) {
      if (data.type == "rec_event2") {
        deathtime <- runif(1, cens.par[1], cens.par[2])
      }
      if (data.type == "rec_event1") {
        deathtime <- rep(runif(ni_, cens.par[1], cens.par[2]), kij_)
      }
      if (data.type == "grouped") {
        deathtime <- runif(sum(kij_), cens.par[1], cens.par[2])
      }
    } else {
      death_rate <- log(2) / cens.par
      if (data.type == "rec_event2") {
        deathtime <- rexp(1, death_rate)
      }
      if (data.type == "rec_event1") {
        deathtime <- rep(rexp(ni_, death_rate), kij_)
      }
      if (data.type == "grouped") {
        deathtime <- runif(sum(kij_), death_rate)
      }
    }

    if (tolower(FUP.type) %in% c("f", "fixed")) {
      cens_time <- pmin(FUP, deathtime)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      cens_time <- pmin(Acc.Dur + FUP - Accrual, deathtime)
    }

    Tijk_H0 <- H0_1(s = -log(U) * exp(-(beta.H0 * Xijk)) / (vi * wij), y = scale_W, p = shape.W)
    Tijk_HA <- H0_1(s = -log(U) * exp(-(beta.HA * Xijk)) / (vi * wij), y = scale_W, p = shape.W)

    ind <- unique(ID_subgroup)
    ind_cum_sub <- sapply(ind, function(i) which(ID_subgroup == i)[length(which(ID_subgroup == i))])

    if ((data.type == "rec_event2") || (data.type == "rec_event1")) {
      if (nb.zerok > 0) {
        Tijk_H0[ind_cum_sub[where.zerok]] <- cens_time[ind_cum_sub[where.zerok]]
        Tijk_HA[ind_cum_sub[where.zerok]] <- cens_time[ind_cum_sub[where.zerok]]
      }
      Yijk_cal_H0 <- pmin(cumsum_m(Tijk_H0, ID_subgroup), cens_time)
      Yijk_cal_HA <- pmin(cumsum_m(Tijk_HA, ID_subgroup), cens_time)

      if (timescale %in% c("calendar", "Calendar") & data.type == "rec_event1") {
        Yijk_H0 <- Yijk_cal_H0
        Yijk_HA <- Yijk_cal_HA
      } else {
        Yijk_H0 <- cumsumrev(Yijk_cal_H0, ID_subgroup)
        Yijk_HA <- cumsumrev(Yijk_cal_HA, ID_subgroup)
      }
      deltaijk <- as.numeric((Yijk_cal_H0 < cens_time))
      deltaijk_ <- as.numeric((Yijk_cal_HA < cens_time))

      Di <- sum(deltaijk)
      Di_ <- sum(deltaijk_)
      Dij <- ave(deltaijk, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      Dij_ <- ave(deltaijk_, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      DiXijk <- sum(deltaijk * Xijk)
      DiXijk_ <- sum(deltaijk_ * Xijk)
      # subs_di <- min(Di+1,ni_)
      # subs_di_ <- min(Di_+1,ni_)
    } else {
      Yijk_H0 <- pmin(Tijk_H0, cens_time)
      Yijk_HA <- pmin(Tijk_HA, cens_time)

      deltaijk <- as.numeric((Yijk_H0 < cens_time))
      deltaijk_ <- as.numeric((Yijk_HA < cens_time))

      Di <- sum(deltaijk)
      Di_ <- sum(deltaijk_)
      Dij <- ave(deltaijk, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      Dij_ <- ave(deltaijk_, ID_subgroup, FUN = cumsum)[ind_cum_sub]
      DiXijk <- sum(deltaijk * Xijk)
      DiXijk_ <- sum(deltaijk_ * Xijk)
    }


    H0Yijk_H0 <- ave(H0(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    H0Yijk_HA <- ave(H0(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    dH0_epsiYijk_H0 <- ave(dH0_epsi(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    dH0_epsiYijk_HA <- ave(dH0_epsi(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    dH0_nuYijk_H0 <- ave(dH0_nu(t = Yijk_H0, y = scale_W, p = shape.W) * exp((beta.H0 * Xijk)), ID_subgroup, FUN = function(x) sum(x, na.rm = T))[ind_cum_sub]
    dH0_nuYijk_HA <- ave(dH0_nu(t = Yijk_HA, y = scale_W, p = shape.W) * exp((beta.HA * Xijk)), ID_subgroup, FUN = function(x) sum(x, na.rm = T))[ind_cum_sub]

    H0Yijk_H0_X <- ave(H0(t = Yijk_H0, y = scale_W, p = shape.W) * Xijk * exp((beta.H0 * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]
    H0Yijk_HA_X <- ave(H0(t = Yijk_HA, y = scale_W, p = shape.W) * Xijk * exp((beta.HA * Xijk)), ID_subgroup, FUN = sum)[ind_cum_sub]

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))


    dbeta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbeta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, H00_X = H0Yijk_H0_X, Dij = Dij)
    }))
    dbeta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbeta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, H00_X = H0Yijk_HA_X, Dij = Dij_)
    }))
    score_beta_H0 <- DiXijk + dbeta_G_THETA_val_H0 / G_THETA_val_H0
    score_beta_HA <- DiXijk_ + dbeta_G_THETA_val_HA / G_THETA_val_HA

    deta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    deta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))

    score_eta_H0 <- digamma(1 / eta) / eta^2 + (log(eta) - 1) / eta^2 + deta_G_THETA_val_H0 / G_THETA_val_H0
    score_eta_HA <- digamma(1 / eta) / eta^2 + (log(eta) - 1) / eta^2 + deta_G_THETA_val_HA / G_THETA_val_HA

    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, Dij = Dij)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, Dij = Dij_)
    }))

    indic <- sum(sapply(c(1:length(kij_)), function(jj) {
      if (Dij[jj] > 0) {
        ln <- c(1:(sum(Dij[jj])))
        return(sum((Dij[jj] - ln) / (1 + theta * (Dij[jj] - ln))))
      } else {
        0
      }
    }))
    indic_ <- sum(sapply(c(1:length(kij_)), function(jj) {
      if (Dij_[jj] > 0) {
        ln_ <- c(1:(sum(Dij_[jj])))
        return(sum((Dij_[jj] - ln_) / (1 + theta * (Dij_[jj] - ln_))))
      } else {
        0
      }
    }))
    score_theta_H0 <- indic + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- indic_ + dtheta_G_THETA_val_HA / G_THETA_val_HA

    depsi_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, dH00 = dH0_epsiYijk_H0, Dij = Dij)
    }))
    depsi_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, dH00 = dH0_epsiYijk_HA, Dij = Dij_)
    }))


    score_epsi_H0 <- -Di * shape.W / scale_W + depsi_G_THETA_val_H0 / G_THETA_val_H0
    score_epsi_HA <- -Di_ * shape.W / scale_W + depsi_G_THETA_val_HA / G_THETA_val_HA


    dnu_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di, eta = eta, theta = theta, H00 = H0Yijk_H0, dH00 = dH0_nuYijk_H0, Dij = Dij)
    }))
    dnu_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0_G_THETA(vi = x, Di = Di_, eta = eta, theta = theta, H00 = H0Yijk_HA, dH00 = dH0_nuYijk_HA, Dij = Dij_)
    }))

    score_nu_H0 <- sum(deltaijk * dh0_nu(Yijk_H0, scale_W, shape.W) / h0(Yijk_H0, scale_W, shape.W), na.rm = TRUE) + dnu_G_THETA_val_H0 / G_THETA_val_H0
    score_nu_HA <- sum(deltaijk_ * dh0_nu(Yijk_HA, scale_W, shape.W) / h0(Yijk_HA, scale_W, shape.W), na.rm = TRUE) + dnu_G_THETA_val_HA / G_THETA_val_HA

    if (npar == 5) {
      score_vector_H0 <- c(score_beta_H0, score_eta_H0, score_theta_H0, score_epsi_H0, score_nu_H0)
      score_vector_HA <- c(score_beta_HA, score_eta_HA, score_theta_HA, score_epsi_HA, score_nu_HA)
    }
    if (sum(is.na(c(score_vector_H0, score_vector_HA))) == 0) {
      U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
      U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)
      sum_score2_H0 <- sum_score2_H0 + U_score2_H0
      sum_score2_HA <- sum_score2_HA + U_score2_HA

      sum_Di_H0 <- sum_Di_H0 + Di / sum(kij_)
      Di_H0 <- Di_H0 + Di
      sum_Di_HA <- sum_Di_HA + Di_ / sum(kij_)
      Di_HA <- Di_HA + Di_
    } else {
      sumNA <- sumNA + 1
    }
  }
  E_score2_H0 <- sum_score2_H0 / (samples.mc - sumNA)
  E_score2_HA <- sum_score2_HA / (samples.mc - sumNA)
  E_Di_H0 <- sum_Di_H0 / (samples.mc - sumNA)
  E_Di_HA <- sum_Di_HA / (samples.mc - sumNA)

  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }

  C_trt <- matrix(data = c(1, 0, 0, 0, 0), nrow = 1, ncol = npar)

  cat("\n")

  rr <- nrow(C_trt)
  if (statistic == "Wald") {
    ncp <- .chi2ncp(alpha * 2, power, df = rr)
    CIC <- t(C_trt %*% (c(beta.HA, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0, 0)))

    Groups <- ceiling(ncp * solve(CIC))

    ncp0 <- (Groups) * t(C_trt %*% (c(beta.H0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(beta.H0, 0, 0, 0, 0)))
    ncpA <- (Groups) * t(C_trt %*% (c(beta.HA, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(beta.HA, 0, 0, 0, 0)))
    estimated.power <- 1 - pchisq(qchisq(1 - alpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA)
  }

  res <- list(
    Groups          = Groups,
    estimated.power = estimated.power,
    events          = c(Groups * Di_H0 / (samples.mc - sumNA), Groups * Di_HA / (samples.mc - sumNA))
    # censoring     = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100), alpha = alpha
  )

  class(res) <- "frailtyDesign"

  return(res)
}






Power_JFM <- function(
    Npts = 400, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
    betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
    shapeR.W = 1, shapeD.W = 1, theta = 0.25, alpha = 1, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  # --- checks ---
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score'")
  }
  if (missing(Npts)) {
    stop("Npts is missing")
  }
  if (Npts <= 0) {
    stop("Npts must be strictly positive")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'unif' or 'pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is max ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & (min(ni) == max(ni))) {
    stop("The upper and lower bounds of ni must be different")
  }
  if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
    stop("The upper and lower bounds of ni must be strictly positive")
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be 'gap' or 'calendar'")
  }
  if ((!betaTest.type %in% c("betaRtest", "betaDtest", "joint"))) {
    stop("betaTest.type must be either 'betaRtest', 'betaDtest' or 'joint' ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }

  # --- setup and parameters ---
  Gc <- Npts
  P_E <- 1 / (1 / ratio + 1)

  theta <- theta
  alpha <- alpha
  betaD.HA <- betaD.HA
  betaD.H0 <- betaD.H0
  betaR.HA <- betaR.HA
  betaR.H0 <- betaR.H0

  npar <- 8
  if (test.type == "1-sided") {
    typeIalpha <- typeIerror
  }
  if (test.type == "2-sided") {
    typeIalpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scaleR_W <- medianR.H0 / log(2)^(1 / shapeR.W)
  scaleD_W <- medianD.H0 / log(2)^(1 / shapeD.W)

  # Helper functions (weibull)
  h0 <- function(t, y, p) {
    (t / y)^(p - 1) * (p / y)
  }
  H0 <- function(t, y, p) {
    (t / y)^p
  }
  dh0_epsi <- function(t, y, p) {
    -p^2 * t^(p - 1) / y^(p + 1)
  }
  dh0_nu <- function(t, y, p) {
    if (length(t) == 0 || all(t == 0)) {
      return(0)
    }
    (1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y))
  }
  dH0_epsi <- function(t, y, p) {
    -p * (t^p) / y^(p + 1)
  }
  dH0_nu <- function(t, y, p) {
    if (length(t) == 0 || all(t == 0)) {
      return(0)
    }
    (t / y)^p * log(t / y)
  }
  H0_1 <- function(s, y, p) {
    y * s^(1 / p)
  }
  # helper func. for gap timescale
  cumsumrev <- function(x) {
    # c( x1, (x2 - x1), (x3 - x2), ...)
    output <- rep(NA, length(x))
    ll <- length(x)
    if (ll > 1) {
      x1 <- x[1:(ll - 1)]
      x2 <- x[2:ll]
      output <- c(x[1], x2 - x1)
    }
    if (ll == 1) {
      output <- x
    }
    return(output)
  }
  # integrand for gamma-frailty part
  G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(omega_i^
      {
        Di + alpha * deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }

  dbetaR_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, H00_R_Z) {
    return(-omega_i^{
      Di + alpha * deltaStar + 1 / theta
    } * H00_R_Z *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dbetaD_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, H00_D_Z) {
    return(-omega_i^{
      Di + alpha * (deltaStar + 1) + 1 / theta - 1
    } * H00_D_Z *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dtheta_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(omega_i^
      {
        Di + alpha * deltaStar + 1 / theta
      } / theta^2 * (1 - log(omega_i) / omega_i) *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dalpha_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(log(omega_i) * omega_i^{
      Di + alpha * deltaStar + 1 / theta - 1
    } * (deltaStar - omega_i^alpha * H00_D) *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dh0R_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, dH00_R) {
    return(-omega_i^{
      Di + alpha * deltaStar + 1 / theta
    } * dH00_R * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dh0D_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, dH00_D) {
    return(-omega_i^{
      Di + alpha * (deltaStar + 1) + 1 / theta - 1
    } * dH00_D * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }

  # --- Accumulators for summations ---
  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- sumNA <- 0
  Di_H0 <- Di_HA <- 0
  sum_deltaStar_H0 <- sum_deltaStar_HA <- 0

  # --- Monte Carlo ---
  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    # sample ni_
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "max"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }

    U_D <- runif(n = 1, 0, 1)

    # random arm allocation, recruitment and individual frailty
    Xi <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
    Accrual <- runif(n = 1, min = 0, max = Acc.Dur)
    omega_i <- rgamma(1, shape = theta^(-1), scale = theta)

    # time to terminal event
    Ti_D_H0 <- H0_1(s = -log(U_D) * exp(-(betaD.H0 * Xi)) / (omega_i^alpha), y = scaleD_W, p = shapeD.W)
    Ti_D_HA <- H0_1(s = -log(U_D) * exp(-(betaD.HA * Xi)) / (omega_i^alpha), y = scaleD_W, p = shapeD.W)

    if (ni_ > 0) {
      U <- runif(n = ni_, 0, 1)
      Tij_H0 <- H0_1(s = -log(U) * exp(-(betaR.H0 * Xi)) / (omega_i), y = scaleR_W, p = shapeR.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(betaR.HA * Xi)) / (omega_i), y = scaleR_W, p = shapeR.W)
    } else {
      # if ni_ = 0 => no rec events
      Tij_H0 <- Ti_D_H0
      Tij_HA <- Ti_D_HA
    }

    # censoring
    if (tolower(FUP.type) %in% c("f", "fixed")) {
      Yi_D_H0 <- pmin(FUP, Ti_D_H0)
      Yi_D_HA <- pmin(FUP, Ti_D_HA)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      Yi_D_H0 <- pmin(Acc.Dur + FUP - Accrual, Ti_D_H0)
      Yi_D_HA <- pmin(Acc.Dur + FUP - Accrual, Ti_D_HA)
    }

    Yij_cal_H0 <- pmin(cumsum(Tij_H0), Yi_D_H0)
    Yij_cal_HA <- pmin(cumsum(Tij_HA), Yi_D_HA)

    deltaij <- as.numeric(Yij_cal_H0 < Yi_D_H0)
    deltaij_ <- as.numeric(Yij_cal_HA < Yi_D_HA)
    Di <- sum(deltaij)
    Di_ <- sum(deltaij_)
    deltaStar <- as.numeric(Ti_D_H0 <= Yi_D_H0)
    deltaStar_ <- as.numeric(Ti_D_HA <= Yi_D_HA)
    DiXi <- sum(deltaij * Xi)
    DiXi_ <- sum(deltaij_ * Xi)

    if (tolower(timescale) == "gap") {
      # Gap scale times
      Yij_H0 <- cumsumrev(Yij_cal_H0)
      Yij_HA <- cumsumrev(Yij_cal_HA)

      # Recurrent hazard/cumhaz in GAP scale = sum of H0(gap_j)
      R0Yij_H0 <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      R0Yij_HA <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      # partial derivatives for shape, scale, etc.
      dR0_epsiYij_H0 <- sum(dH0_epsi(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      dR0_epsiYij_HA <- sum(dH0_epsi(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_nuYij_H0 <- sum(dH0_nu(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.H0 * Xi)
      dR0_nuYij_HA <- sum(dH0_nu(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.HA * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      R0Yij_H0_X <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * Xi * exp(betaR.H0 * Xi)
      R0Yij_HA_X <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * Xi * exp(betaR.HA * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    } else {
      Yij_H0 <- Yij_cal_H0
      Yij_HA <- Yij_cal_HA
      ## ----------------------------------------------------
      ## CALENDAR TIMESCALE, using telescoping:
      ##    sum_j [ R0(T_{ij}) - R0(T_{i,j-1}) ] = R0(T_{i,n_i}) - R0(Ti,0).
      ##    Note: Ti_ni = T_max_H0
      ##          Ti,0  = T_min_H0 (= 0 if there's no left truncated data).
      ## ----------------------------------------------------

      # Number of recurrent events (calendar time)
      n_j_H0 <- length(Yij_cal_H0) # Yij_cal_H0 = sorted event/censor times for H0 scenario
      n_j_HA <- length(Yij_cal_HA) # Yij_cal_HA = likewise for HA scenario

      ## ------- H0 scenario -------
      if (n_j_H0 > 0) {
        # T_max_H0 = largest calendar event time
        T_max_H0 <- Yij_cal_H0[n_j_H0]
        T_min_H0 <- Yij_cal_H0[1]

        # Total 'cumHaz' from 0 to T_max_H0:
        R0Yij_H0 <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        # Derivatives wrt scale (epsi) and shape (nu):
        dR0_epsiYij_H0 <- (
          dH0_epsi(T_max_H0, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_nuYij_H0 <- (
          dH0_nu(T_max_H0, scaleR_W, shapeR.W) -
            dH0_nu(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        # If Xi>0, partial wrt betaR "pulls down" factor Xi:
        R0Yij_H0_X <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.H0 * Xi)
      } else {
        # No recurrent events => 0 contribution
        R0Yij_H0 <- 0
        dR0_epsiYij_H0 <- 0
        dR0_nuYij_H0 <- 0
        R0Yij_H0_X <- 0
      }

      # Terminal event part (Weibull) remains unchanged:
      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)

      ## ------- HA scenario -------
      if (n_j_HA > 0) {
        T_max_HA <- Yij_cal_HA[n_j_HA]
        T_min_HA <- Yij_cal_HA[1]
        R0Yij_HA <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_epsiYij_HA <- (
          dH0_epsi(T_max_HA, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_nuYij_HA <- (
          dH0_nu(T_max_HA, scaleR_W, shapeR.W) -
            dH0_nu(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        R0Yij_HA_X <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.HA * Xi)
      } else {
        R0Yij_HA <- 0
        dR0_epsiYij_HA <- 0
        dR0_nuYij_HA <- 0
        R0Yij_HA_X <- 0
      }

      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    } ## end of 'else' calendar timescale

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA)
    }))

    # partial wrt betaR
    dbetaR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, R0Yij_H0_X)
    }))
    dbetaR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, R0Yij_HA_X)
    }))
    score_betaR_H0 <- DiXi + dbetaR_G_THETA_val_H0 / G_THETA_val_H0
    score_betaR_HA <- DiXi_ + dbetaR_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt betaD
    dbetaD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, H0Yi_H0_X)
    }))
    dbetaD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, H0Yi_HA_X)
    }))
    score_betaD_H0 <- (deltaStar * Xi) + dbetaD_G_THETA_val_H0 / G_THETA_val_H0
    score_betaD_HA <- (deltaStar_ * Xi) + dbetaD_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt theta
    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA)
    }))
    score_theta_H0 <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt alpha
    dalpha_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dalpha_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0)
    }))
    dalpha_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dalpha_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA)
    }))
    score_alpha_H0 <- dalpha_G_THETA_val_H0 / G_THETA_val_H0
    score_alpha_HA <- dalpha_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt scaleR => epsiR (Weibull)
    depsiR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, dR0_epsiYij_H0)
    }))
    depsiR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, dR0_epsiYij_HA)
    }))
    score_epsiR_H0 <- -Di * shapeR.W / scaleR_W + depsiR_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiR_HA <- -Di_ * shapeR.W / scaleR_W + depsiR_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt nuR => shape param for rec
    dnuR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, dR0_nuYij_H0)
    }))
    dnuR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, dR0_nuYij_HA)
    }))
    # # the 'observed' part for shape is sum_j [deltaij * (dh0_nu(...) / h0(...))]
    score_nuR_H0 <- sum(deltaij * dh0_nu(Yij_H0, scaleR_W, shapeR.W) / h0(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) +
      dnuR_G_THETA_val_H0 / G_THETA_val_H0
    score_nuR_HA <- sum(deltaij_ * dh0_nu(Yij_HA, scaleR_W, shapeR.W) / h0(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) +
      dnuR_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt scaleD => epsiD
    depsiD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, dH00_D = dH0_epsiYi_H0)
    }))
    depsiD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, dH00_D = dH0_epsiYi_HA)
    }))
    score_epsiD_H0 <- -deltaStar * shapeD.W / scaleD_W + depsiD_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiD_HA <- -deltaStar_ * shapeD.W / scaleD_W + depsiD_G_THETA_val_HA / G_THETA_val_HA

    # partial wrt nuD => shape param for death
    dnuD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(x, Di, deltaStar, theta, alpha, R0Yij_H0, H0Yi_H0, dH00_D = dH0_nuYi_H0)
    }))
    dnuD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(x, Di_, deltaStar_, theta, alpha, R0Yij_HA, H0Yi_HA, dH00_D = dH0_nuYi_HA)
    }))
    score_nuD_H0 <- deltaStar * dh0_nu(Yi_D_H0, scaleD_W, shapeD.W) / h0(Yi_D_H0, scaleD_W, shapeD.W) +
      dnuD_G_THETA_val_H0 / G_THETA_val_H0
    score_nuD_HA <- deltaStar_ * dh0_nu(Yi_D_HA, scaleD_W, shapeD.W) / h0(Yi_D_HA, scaleD_W, shapeD.W) +
      dnuD_G_THETA_val_HA / G_THETA_val_HA

    # Build full score vector
    score_vector_H0 <- c(
      score_betaR_H0, score_betaD_H0, score_theta_H0,
      score_alpha_H0, score_epsiR_H0, score_nuR_H0,
      score_epsiD_H0, score_nuD_H0
    )
    score_vector_HA <- c(
      score_betaR_HA, score_betaD_HA, score_theta_HA,
      score_alpha_HA, score_epsiR_HA, score_nuR_HA,
      score_epsiD_HA, score_nuD_HA
    )

    # update sums if no NA
    if (!anyNA(c(score_vector_H0, score_vector_HA))) {
      U_score2_H0 <- score_vector_H0 %*% t(score_vector_H0)
      U_score2_HA <- score_vector_HA %*% t(score_vector_HA)

      sum_score2_H0 <- sum_score2_H0 + U_score2_H0
      sum_score2_HA <- sum_score2_HA + U_score2_HA

      if (ni_ > 0) {
        sum_Di_H0 <- sum_Di_H0 + Di / ni_
        Di_H0 <- Di_H0 + Di
        sum_Di_HA <- sum_Di_HA + Di_ / ni_
        Di_HA <- Di_HA + Di_
      }
      sum_deltaStar_H0 <- sum_deltaStar_H0 + deltaStar
      sum_deltaStar_HA <- sum_deltaStar_HA + deltaStar_
    } else {
      sumNA <- sumNA + 1
    }
  } # end of MC


  # --- Averages & inverse info ---
  n_eff <- (samples.mc - sumNA)
  E_score2_H0 <- sum_score2_H0 / n_eff
  E_score2_HA <- sum_score2_HA / n_eff
  E_Di_H0 <- sum_Di_H0 / n_eff
  E_Di_HA <- sum_Di_HA / n_eff
  E_deltaStar_H0 <- sum_deltaStar_H0 / n_eff
  E_deltaStar_HA <- sum_deltaStar_HA / n_eff

  if (is.square.matrix(E_score2_H0) && is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) && is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }

  C_trt <- switch(betaTest.type,
    "betaRtest" = matrix(c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1),
    "betaDtest" = matrix(c(0, 1, 0, 0, 0, 0, 0, 0), nrow = 1),
    "joint" = matrix(
      c(
        1, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0
      ),
      nrow = 2, byrow = TRUE
    )
  )

  if (statistic == "Wald") {
    # Evaluate the non-central parameters (=ncp)
    ref0 <- c(betaR.H0, betaD.H0, 0, 0, 0, 0, 0, 0)
    refA <- c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)

    ncp0 <- Gc * t(C_trt %*% ref0) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% ref0)
    ncpA <- Gc * t(C_trt %*% refA) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% refA)

    rr <- nrow(C_trt) # nb. constraints
    crit_val <- qchisq(1 - typeIalpha * 2, df = rr, ncp = ncp0)
    pow_val <- 1 - pchisq(crit_val, df = rr, ncp = ncpA)
  } else {
    stop("Score test not implemented.")
  }

  res <- list(
    estimated.power         = pow_val,
    events.rec    = c(Gc * Di_H0 / samples.mc, Gc * Di_HA / samples.mc),
    events.D      = c(Gc * sum_deltaStar_H0 / samples.mc, Gc * sum_deltaStar_HA / samples.mc)
    # censoring.rec = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100),
    # censoring.D = c((1 - E_deltaStar_H0) * 100, (1 - E_deltaStar_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}

Nsn_JFM <- function(
    power = 0.8, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
    betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
    shapeR.W = 1, shapeD.W = 1, theta = 0.25, alpha = 1, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  if (missing(power)) {
    stop("Power is missing")
  }
  if ((power <= 0) | (power > 0.99)) {
    stop("power must be between ]0,1[")
  }
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'Unif' or 'Pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is max ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & (min(ni) == max(ni))) {
    stop("The upper and lower bounds of ni must be different")
  }
  if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
    stop("The upper and lower bounds of ni must be strictly positive")
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be 'gap' or 'calendar'")
  }
  if ((!betaTest.type %in% c("betaRtest", "betaDtest", "joint"))) {
    stop("betaTest.type must be either 'betaRtest', 'betaDtest' or 'joint' ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }


  P_E <- 1 / (1 / ratio + 1)

  theta <- theta
  alpha <- alpha
  betaD.HA <- betaD.HA
  betaD.H0 <- betaD.H0
  betaR.HA <- betaR.HA
  betaR.H0 <- betaR.H0

  npar <- 8
  if (test.type == "1-sided") {
    typeIalpha <- typeIerror
  }
  if (test.type == "2-sided") {
    typeIalpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scaleR_W <- medianR.H0 / log(2)^(1 / shapeR.W)
  scaleD_W <- medianD.H0 / log(2)^(1 / shapeD.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }
  cumsumrev <- function(x) {
    if (length(x) == 0) {
      return(x)
    }
    c(x[1], diff(x))
  }
  G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(omega_i^
      {
        Di + alpha * deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dbetaR_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, H00_R_Z) {
    return(-omega_i^{
      Di + alpha * deltaStar + 1 / theta
    } * H00_R_Z *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dbetaD_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, H00_D_Z) {
    return(-omega_i^{
      Di + alpha * (deltaStar + 1) + 1 / theta - 1
    } * H00_D_Z *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dtheta_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(omega_i^
      {
        Di + alpha * deltaStar + 1 / theta
      } / theta^2 * (1 - log(omega_i) / omega_i) *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dalpha_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D) {
    return(log(omega_i) * omega_i^{
      Di + alpha * deltaStar + 1 / theta - 1
    } * (deltaStar - omega_i^alpha * H00_D) *
      exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dh0R_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, dH00_R) {
    return(-omega_i^{
      Di + alpha * deltaStar + 1 / theta
    } * dH00_R * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }
  dh0D_G_THETA <- function(omega_i, Di, deltaStar, theta, alpha, H00_R, H00_D, dH00_D) {
    return(-omega_i^{
      Di + alpha * (deltaStar + 1) + 1 / theta - 1
    } * dH00_D * exp(-omega_i * (1 / theta - 1 + H00_R) - omega_i^alpha * H00_D))
  }

  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- sumNA <- 0
  Di_H0 <- Di_HA <- 0
  sum_deltaStar_H0 <- sum_deltaStar_HA <- 0

  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "max"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }

    U_D <- runif(n = 1, 0, 1)

    Xi <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
    Accrual <- runif(n = 1, min = 0, max = Acc.Dur)
    omega_i <- rgamma(1, shape = theta^(-1), scale = theta)

    Ti_D_H0 <- H0_1(s = -log(U_D) * exp(-(betaD.H0 * Xi)) / (omega_i^alpha), y = scaleD_W, p = shapeD.W)
    Ti_D_HA <- H0_1(s = -log(U_D) * exp(-(betaD.HA * Xi)) / (omega_i^alpha), y = scaleD_W, p = shapeD.W)

    if (ni_ > 0) {
      U <- runif(n = ni_, 0, 1)
      Tij_H0 <- H0_1(s = -log(U) * exp(-(betaR.H0 * Xi)) / (omega_i), y = scaleR_W, p = shapeR.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(betaR.HA * Xi)) / (omega_i), y = scaleR_W, p = shapeR.W)
    } else {
      Tij_H0 <- Ti_D_H0
      Tij_HA <- Ti_D_HA
    }

    if (tolower(FUP.type) %in% c("f", "fixed")) {
      Yi_D_H0 <- pmin(FUP, Ti_D_H0)
      Yi_D_HA <- pmin(FUP, Ti_D_HA)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      Yi_D_H0 <- pmin(Acc.Dur + FUP - Accrual, Ti_D_H0)
      Yi_D_HA <- pmin(Acc.Dur + FUP - Accrual, Ti_D_HA)
    }

    Yij_cal_H0 <- pmin(cumsum(Tij_H0), Yi_D_H0)
    Yij_cal_HA <- pmin(cumsum(Tij_HA), Yi_D_HA)

    deltaij <- as.numeric(Yij_cal_H0 < Yi_D_H0)
    deltaij_ <- as.numeric(Yij_cal_HA < Yi_D_HA)
    Di <- sum(deltaij)
    Di_ <- sum(deltaij_)
    deltaStar <- as.numeric(Ti_D_H0 <= Yi_D_H0)
    deltaStar_ <- as.numeric(Ti_D_HA <= Yi_D_HA)
    DiXi <- sum(deltaij * Xi)
    DiXi_ <- sum(deltaij_ * Xi)

    if (tolower(timescale) == "gap") {
      Yij_H0 <- cumsumrev(Yij_cal_H0)
      Yij_HA <- cumsumrev(Yij_cal_HA)

      R0Yij_H0 <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      R0Yij_HA <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_epsiYij_H0 <- sum(dH0_epsi(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      dR0_epsiYij_HA <- sum(dH0_epsi(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_nuYij_H0 <- sum(dH0_nu(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.H0 * Xi)
      dR0_nuYij_HA <- sum(dH0_nu(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.HA * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      R0Yij_H0_X <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * Xi * exp(betaR.H0 * Xi)
      R0Yij_HA_X <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * Xi * exp(betaR.HA * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    } else {
      Yij_H0 <- Yij_cal_H0
      Yij_HA <- Yij_cal_HA
      n_j_H0 <- length(Yij_cal_H0)
      n_j_HA <- length(Yij_cal_HA)

      if (n_j_H0 > 0) {
        T_max_H0 <- Yij_cal_H0[n_j_H0]
        T_min_H0 <- Yij_cal_H0[1]

        R0Yij_H0 <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_epsiYij_H0 <- (
          dH0_epsi(T_max_H0, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_nuYij_H0 <- (
          dH0_nu(T_max_H0, scaleR_W, shapeR.W) -
            dH0_nu(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        R0Yij_H0_X <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.H0 * Xi)
      } else {
        R0Yij_H0 <- 0
        dR0_epsiYij_H0 <- 0
        dR0_nuYij_H0 <- 0
        R0Yij_H0_X <- 0
      }

      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)

      if (n_j_HA > 0) {
        T_max_HA <- Yij_cal_HA[n_j_HA]
        T_min_HA <- Yij_cal_HA[1]
        R0Yij_HA <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_epsiYij_HA <- (
          dH0_epsi(T_max_HA, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_nuYij_HA <- (
          dH0_nu(T_max_HA, scaleR_W, shapeR.W) -
            dH0_nu(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        R0Yij_HA_X <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.HA * Xi)
      } else {
        R0Yij_HA <- 0
        dR0_epsiYij_HA <- 0
        dR0_nuYij_HA <- 0
        R0Yij_HA_X <- 0
      }

      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    }

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    dbetaR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_R_Z = R0Yij_H0_X)
    }))
    dbetaR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_R_Z = R0Yij_HA_X)
    }))
    score_betaR_H0 <- DiXi + dbetaR_G_THETA_val_H0 / G_THETA_val_H0
    score_betaR_HA <- DiXi_ + dbetaR_G_THETA_val_HA / G_THETA_val_HA

    dbetaD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_D_Z = H0Yi_H0_X)
    }))
    dbetaD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_D_Z = H0Yi_HA_X)
    }))
    score_betaD_H0 <- deltaStar * Xi + dbetaD_G_THETA_val_H0 / G_THETA_val_H0
    score_betaD_HA <- deltaStar_ * Xi + dbetaD_G_THETA_val_HA / G_THETA_val_HA

    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    score_theta_H0 <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_HA / G_THETA_val_HA

    dalpha_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dalpha_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    dalpha_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dalpha_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    score_alpha_H0 <- dalpha_G_THETA_val_H0 / G_THETA_val_H0
    score_alpha_HA <- dalpha_G_THETA_val_HA / G_THETA_val_HA

    depsiR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_epsiYij_H0)
    }))
    depsiR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_epsiYij_HA)
    }))

    score_epsiR_H0 <- -Di * shapeR.W / scaleR_W + depsiR_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiR_HA <- -Di_ * shapeR.W / scaleR_W + depsiR_G_THETA_val_HA / G_THETA_val_HA


    dnuR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_nuYij_H0)
    }))
    dnuR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_nuYij_HA)
    }))

    score_nuR_H0 <- sum(deltaij * dh0_nu(Yij_H0, scaleR_W, shapeR.W) / h0(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_H0 / G_THETA_val_H0
    score_nuR_HA <- sum(deltaij_ * dh0_nu(Yij_HA, scaleR_W, shapeR.W) / h0(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_HA / G_THETA_val_HA

    depsiD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_epsiYi_H0)
    }))
    depsiD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_epsiYi_HA)
    }))


    score_epsiD_H0 <- -deltaStar * shapeD.W / scaleD_W + depsiD_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiD_HA <- -deltaStar_ * shapeD.W / scaleD_W + depsiD_G_THETA_val_HA / G_THETA_val_HA

    dnuD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, alpha = alpha, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_nuYi_H0)
    }))
    dnuD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, alpha = alpha, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_nuYi_HA)
    }))

    score_nuD_H0 <- deltaStar * dh0_nu(Yi_D_H0, scaleD_W, shapeD.W) / h0(Yi_D_H0, scaleD_W, shapeD.W) + dnuD_G_THETA_val_H0 / G_THETA_val_H0
    score_nuD_HA <- deltaStar_ * dh0_nu(Yi_D_HA, scaleD_W, shapeD.W) / h0(Yi_D_HA, scaleD_W, shapeD.W) + dnuD_G_THETA_val_HA / G_THETA_val_HA

    if (npar == 8) {
      score_vector_H0 <- c(score_betaR_H0, score_betaD_H0, score_theta_H0, score_alpha_H0, score_epsiR_H0, score_nuR_H0, score_epsiD_H0, score_nuD_H0)
      score_vector_HA <- c(score_betaR_HA, score_betaD_HA, score_theta_HA, score_alpha_HA, score_epsiR_HA, score_nuR_HA, score_epsiD_HA, score_nuD_HA)
    }
    if (sum(is.na(c(score_vector_H0, score_vector_HA))) == 0) {
      U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
      U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)

      sum_score2_H0 <- sum_score2_H0 + U_score2_H0
      sum_score2_HA <- sum_score2_HA + U_score2_HA

      if (ni_ > 0) {
        sum_Di_H0 <- sum_Di_H0 + Di / ni_
        Di_H0 <- Di_H0 + Di
        sum_Di_HA <- sum_Di_HA + Di_ / ni_
        Di_HA <- Di_HA + Di_
      }
      sum_deltaStar_H0 <- sum_deltaStar_H0 + deltaStar
      sum_deltaStar_HA <- sum_deltaStar_HA + deltaStar_
    } else {
      sumNA <- sumNA + 1
    }
  }
  E_score2_H0 <- sum_score2_H0 / (samples.mc - sumNA)
  E_score2_HA <- sum_score2_HA / (samples.mc - sumNA)
  E_Di_H0 <- sum_Di_H0 / (samples.mc - sumNA)
  E_Di_HA <- sum_Di_HA / (samples.mc - sumNA)
  E_deltaStar_H0 <- sum_deltaStar_H0 / (samples.mc - sumNA)
  E_deltaStar_HA <- sum_deltaStar_HA / (samples.mc - sumNA)

  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }

  C_trt <- switch(betaTest.type,
    "betaRtest" = matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "betaDtest" = matrix(data = c(0, 1, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "joint" = matrix(data = c(
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 2, ncol = npar)
  )

  cat("\n")
  rr <- nrow(C_trt)

  if (statistic == "Wald") {
    ncp <- .chi2ncp(typeIalpha * 2, power, df = rr)
    CIC <- t(C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)))

    Npts <- ceiling(ncp * solve(CIC))

    ncp0 <- (Npts) * t(C_trt %*% (c(betaR.H0, betaD.H0, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(betaR.H0, betaD.H0, 0, 0, 0, 0, 0, 0)))
    ncpA <- (Npts) * t(C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)))
    estimated.power <- 1 - pchisq(qchisq(1 - typeIalpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA)
  }

  res <- list(
    model = "JFM",
    method = "ssize",
    timescale = timescale,
    Npts = Npts,
    estimated.power = estimated.power,
    events.rec = c(Npts * Di_H0 / samples.mc, Npts * Di_HA / samples.mc),
    events.D = c(Npts * sum_deltaStar_H0 / samples.mc, Npts * sum_deltaStar_HA / samples.mc)
    # censoring.rec = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100),
    # censoring.D = c((1 - E_deltaStar_H0) * 100, (1 - E_deltaStar_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}







Power_GJFM <- function(
    Npts = 400, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
    betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
    shapeR.W = 1, shapeD.W = 1, theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(Npts)) {
    stop("Npts is missing")
  }
  if (Npts <= 0) {
    stop("Npts must be strictly positive")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'Unif' or 'Pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is max ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & (min(ni) == max(ni))) {
    stop("The upper and lower bounds of ni must be different")
  }
  if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
    stop("The upper and lower bounds of ni must be strictly positive")
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be 'gap' or 'calendar'")
  }
  if ((!betaTest.type %in% c("betaRtest", "betaDtest", "joint"))) {
    stop("betaTest.type must be either 'betaRtest', 'betaDtest' or 'joint' ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }


  Gc <- Npts
  P_E <- 1 / (1 / ratio + 1)

  theta <- theta
  eta <- eta
  betaD.HA <- betaD.HA
  betaD.H0 <- betaD.H0
  betaR.HA <- betaR.HA
  betaR.H0 <- betaR.H0

  npar <- 8
  if (test.type == "1-sided") {
    typeIalpha <- typeIerror
  }
  if (test.type == "2-sided") {
    typeIalpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scaleR_W <- medianR.H0 / log(2)^(1 / shapeR.W)
  scaleD_W <- medianD.H0 / log(2)^(1 / shapeD.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }
  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }
  cumsumrev <- function(x) {
    output <- rep(NA, length(x))
    ll <- length(x)
    if (ll > 1) {
      x1 <- x[1:(ll - 1)]
      x2 <- x[2:ll]
      output <- c(x[1], x2 - x1)
    }
    if (ll == 1) {
      output <- x
    }
    return(output)
  }
  G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    return(omega_i^
      {
        Di + deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }

  dbetaR_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, H00_R_Z) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * ((Di * eta + 1) * H00_R_Z / (1 + eta * omega_i * H00_R)^(Di + 1 / eta + 1)))
  }

  dbetaD_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, H00_D_Z) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * H00_D_Z * exp(-omega_i * (1 / theta - 1 + H00_D)) /
      (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }


  dtheta_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    return(omega_i^
      {
        Di + deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * omega_i * H00_R)^(Di + 1 / eta) *
      (omega_i - log(omega_i)) / theta^2)
  }

  deta_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    alpha <- omega_i * H00_R

    deriveta <- exp((Di + 1 / eta) * log(1 + eta * alpha)) * ((-1 / eta^2) * log(1 + eta * alpha) + (Di + 1 / eta) * alpha / (1 + eta * alpha))
    return(-omega_i^{
      Di + deltaStar + 1 / theta - 1
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * alpha)^{
      2 * (Di + 1 / eta)
    } *
      deriveta)
  }


  dh0R_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, dH00_R) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta - 1
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * (Di + 1 / eta) * eta * omega_i * dH00_R / (1 + eta * omega_i * H00_R)^(Di + 1 / eta + 1))
  }

  dh0D_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, dH00_D) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * dH00_D / (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }

  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- 0
  Di_H0 <- Di_HA <- 0
  sum_deltaStar_H0 <- sum_deltaStar_HA <- 0

  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "max"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }

    U_D <- runif(n = 1, 0, 1)

    Xi <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
    Accrual <- runif(n = 1, min = 0, max = Acc.Dur)

    omega_i <- rgamma(1, shape = theta^(-1), scale = theta)
    vi <- rgamma(1, shape = eta^(-1), scale = eta)


    Ti_D_H0 <- H0_1(s = -log(U_D) * exp(-(betaD.H0 * Xi)) / (omega_i), y = scaleD_W, p = shapeD.W)
    Ti_D_HA <- H0_1(s = -log(U_D) * exp(-(betaD.HA * Xi)) / (omega_i), y = scaleD_W, p = shapeD.W)


    if (ni_ > 0) {
      U <- runif(n = ni_, 0, 1)
      Tij_H0 <- H0_1(s = -log(U) * exp(-(betaR.H0 * Xi)) / (omega_i * vi), y = scaleR_W, p = shapeR.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(betaR.HA * Xi)) / (omega_i * vi), y = scaleR_W, p = shapeR.W)
    } else {
      Tij_H0 <- Ti_D_H0
      Tij_HA <- Ti_D_HA
    }

    if (tolower(FUP.type) %in% c("f", "fixed")) {
      Yi_D_H0 <- pmin(FUP, Ti_D_H0)
      Yi_D_HA <- pmin(FUP, Ti_D_HA)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      Yi_D_H0 <- pmin(Acc.Dur + FUP - Accrual, Ti_D_H0)
      Yi_D_HA <- pmin(Acc.Dur + FUP - Accrual, Ti_D_HA)
    }

    Yij_cal_H0 <- pmin(cumsum(Tij_H0), Yi_D_H0)
    Yij_cal_HA <- pmin(cumsum(Tij_HA), Yi_D_HA)

    deltaij <- as.numeric(Yij_cal_H0 < Yi_D_H0)
    deltaij_ <- as.numeric(Yij_cal_HA < Yi_D_HA)
    Di <- sum(deltaij)
    Di_ <- sum(deltaij_)
    deltaStar <- as.numeric(Ti_D_H0 <= Yi_D_H0)
    deltaStar_ <- as.numeric(Ti_D_HA <= Yi_D_HA)
    DiXi <- sum(deltaij * Xi)
    DiXi_ <- sum(deltaij_ * Xi)

    if (tolower(timescale) == "gap") {
      Yij_H0 <- cumsumrev(Yij_cal_H0)
      Yij_HA <- cumsumrev(Yij_cal_HA)

      R0Yij_H0 <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      R0Yij_HA <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_epsiYij_H0 <- sum(dH0_epsi(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      dR0_epsiYij_HA <- sum(dH0_epsi(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_nuYij_H0 <- sum(dH0_nu(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.H0 * Xi)
      dR0_nuYij_HA <- sum(dH0_nu(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.HA * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      R0Yij_H0_X <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * Xi * exp(betaR.H0 * Xi)
      R0Yij_HA_X <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * Xi * exp(betaR.HA * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    } else {
      Yij_H0 <- Yij_cal_H0
      Yij_HA <- Yij_cal_HA
      n_j_H0 <- length(Yij_cal_H0)
      n_j_HA <- length(Yij_cal_HA)

      if (n_j_H0 > 0) {
        T_max_H0 <- Yij_cal_H0[n_j_H0]
        T_min_H0 <- Yij_cal_H0[1]

        R0Yij_H0 <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_epsiYij_H0 <- (
          dH0_epsi(T_max_H0, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_nuYij_H0 <- (
          dH0_nu(T_max_H0, scaleR_W, shapeR.W) -
            dH0_nu(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        R0Yij_H0_X <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.H0 * Xi)
      } else {
        R0Yij_H0 <- 0
        dR0_epsiYij_H0 <- 0
        dR0_nuYij_H0 <- 0
        R0Yij_H0_X <- 0
      }

      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)

      if (n_j_HA > 0) {
        T_max_HA <- Yij_cal_HA[n_j_HA]
        T_min_HA <- Yij_cal_HA[1]
        R0Yij_HA <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_epsiYij_HA <- (
          dH0_epsi(T_max_HA, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_nuYij_HA <- (
          dH0_nu(T_max_HA, scaleR_W, shapeR.W) -
            dH0_nu(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        R0Yij_HA_X <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.HA * Xi)
      } else {
        R0Yij_HA <- 0
        dR0_epsiYij_HA <- 0
        dR0_nuYij_HA <- 0
        R0Yij_HA_X <- 0
      }

      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    }

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    dbetaR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_R_Z = R0Yij_H0_X)
    }))
    dbetaR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_R_Z = R0Yij_HA_X)
    }))
    score_betaR_H0 <- DiXi + dbetaR_G_THETA_val_H0 / G_THETA_val_H0
    score_betaR_HA <- DiXi_ + dbetaR_G_THETA_val_HA / G_THETA_val_HA

    dbetaD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_D_Z = H0Yi_H0_X)
    }))
    dbetaD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_D_Z = H0Yi_HA_X)
    }))
    score_betaD_H0 <- deltaStar * Xi + dbetaD_G_THETA_val_H0 / G_THETA_val_H0
    score_betaD_HA <- deltaStar_ * Xi + dbetaD_G_THETA_val_HA / G_THETA_val_HA

    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    score_theta_H0 <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_HA / G_THETA_val_HA

    deta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    deta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))
    l <- l_ <- 0
    if (Di > 0) {
      l <- 1:Di
    }
    if (Di_ > 0) {
      l_ <- 1:Di_
    }

    score_eta_H0 <- Di / eta - sum((1 / eta^2) / (Di + 1 / eta - l)) + deta_G_THETA_val_H0 / G_THETA_val_H0
    score_eta_HA <- Di_ / eta - sum((1 / eta^2) / (Di_ + 1 / eta - l_)) + deta_G_THETA_val_HA / G_THETA_val_HA

    depsiR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_epsiYij_H0)
    }))
    depsiR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_epsiYij_HA)
    }))

    score_epsiR_H0 <- -Di * shapeR.W / scaleR_W + depsiR_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiR_HA <- -Di_ * shapeR.W / scaleR_W + depsiR_G_THETA_val_HA / G_THETA_val_HA

    dnuR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_nuYij_H0)
    }))
    dnuR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_nuYij_HA)
    }))

    score_nuR_H0 <- sum(deltaij * dh0_nu(Yij_H0, scaleR_W, shapeR.W) / h0(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_H0 / G_THETA_val_H0
    score_nuR_HA <- sum(deltaij_ * dh0_nu(Yij_HA, scaleR_W, shapeR.W) / h0(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_HA / G_THETA_val_HA


    depsiD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_epsiYi_H0)
    }))
    depsiD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_epsiYi_HA)
    }))

    score_epsiD_H0 <- -deltaStar * shapeD.W / scaleD_W + depsiD_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiD_HA <- -deltaStar_ * shapeD.W / scaleD_W + depsiD_G_THETA_val_HA / G_THETA_val_HA

    dnuD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_nuYi_H0)
    }))
    dnuD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_nuYi_HA)
    }))

    score_nuD_H0 <- deltaStar * dh0_nu(Yi_D_H0, scaleD_W, shapeD.W) / h0(Yi_D_H0, scaleD_W, shapeD.W) + dnuD_G_THETA_val_H0 / G_THETA_val_H0
    score_nuD_HA <- deltaStar_ * dh0_nu(Yi_D_HA, scaleD_W, shapeD.W) / h0(Yi_D_HA, scaleD_W, shapeD.W) + dnuD_G_THETA_val_HA / G_THETA_val_HA

    if (npar == 8) {
      score_vector_H0 <- c(score_betaR_H0, score_betaD_H0, score_theta_H0, score_eta_H0, score_epsiR_H0, score_nuR_H0, score_epsiD_H0, score_nuD_H0)
      score_vector_HA <- c(score_betaR_HA, score_betaD_HA, score_theta_HA, score_eta_HA, score_epsiR_HA, score_nuR_HA, score_epsiD_HA, score_nuD_HA)
    }

    U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
    U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)

    sum_score2_H0 <- sum_score2_H0 + U_score2_H0
    sum_score2_HA <- sum_score2_HA + U_score2_HA

    if (ni_ > 0) {
      sum_Di_H0 <- sum_Di_H0 + Di / ni_
      Di_H0 <- Di_H0 + Di
      sum_Di_HA <- sum_Di_HA + Di_ / ni_
      Di_HA <- Di_HA + Di_
    }
    sum_deltaStar_H0 <- sum_deltaStar_H0 + deltaStar
    sum_deltaStar_HA <- sum_deltaStar_HA + deltaStar_
  }
  E_score2_H0 <- sum_score2_H0 / samples.mc
  E_score2_HA <- sum_score2_HA / samples.mc
  E_Di_H0 <- sum_Di_H0 / samples.mc
  E_Di_HA <- sum_Di_HA / samples.mc
  E_deltaStar_H0 <- sum_deltaStar_H0 / samples.mc
  E_deltaStar_HA <- sum_deltaStar_HA / samples.mc
  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }


  C_trt <- switch(betaTest.type,
    "betaRtest" = matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "betaDtest" = matrix(data = c(0, 1, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "joint" = matrix(data = c(
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 2, ncol = npar)
  )

  rr <- nrow(C_trt)

  cat("\n")

  if (statistic == "Wald") {
    ncp0 <- (Gc) * t(C_trt %*% (c(betaR.H0, betaD.H0, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_H0 %*% t(C_trt)) %*% (C_trt %*% (c(betaR.H0, betaD.H0, 0, 0, 0, 0, 0, 0)))
    ncpA <- (Gc) * t(C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)))
  }

  res <- list(
    estimated.power         = 1 - pchisq(qchisq(1 - typeIalpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA),
    events.rec    = c(Gc * Di_H0 / samples.mc, Gc * Di_HA / samples.mc),
    events.D      = c(Gc * sum_deltaStar_H0 / samples.mc, Gc * sum_deltaStar_HA / samples.mc)
    # censoring.rec = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100),
    # censoring.D = c((1 - E_deltaStar_H0) * 100, (1 - E_deltaStar_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}

Nsn_GJFM <- function(
    power = 0.8, ni = 8, ni.type = "max",
    Acc.Dur = 0, FUP = 12, FUP.type = "UpToEnd", medianR.H0 = 3, medianD.H0 = 10,
    betaTest.type = "joint", betaR.H0 = 0, betaR.HA = log(0.75), betaD.H0 = 0, betaD.HA = log(0.85),
    shapeR.W = 1, shapeD.W = 1, theta = 0.25, eta = 0.5, ratio = 1, samples.mc = 1e4, seed = 42,
    timescale = "gap", statistic = "Wald", typeIerror = 0.05, test.type = "2-sided") {
  if ((!test.type %in% c("1-sided", "2-sided"))) {
    stop("test.type must be '1-sided' or '2-sided'")
  }
  if ((!statistic %in% c("Wald", "Score"))) {
    stop("statistic must be either 'Wald' or 'Score' ")
  }
  if (missing(power)) {
    stop("power is missing")
  }
  if (power <= 0) {
    stop("Npts must be strictly positive")
  }
  if ((typeIerror <= 0) | (typeIerror > 0.2)) {
    stop("We recommend a typeIerror in the interval (0, 0.2)")
  }
  if (theta <= 0) {
    stop("theta must be strictly positive")
  }
  if (eta <= 0) {
    stop("eta must be strictly positive")
  }
  if ((!tolower(FUP.type) %in% c("fixed", "uptoend", "end"))) {
    stop("FUP.type must be either 'Fixed' or 'UptoEnd' ")
  }
  if ((!tolower(ni.type) %in% c("m", "max", "maximum", "u", "uniform", "unif", "p", "pois", "poisson"))) {
    stop("ni.type must be either 'max' or 'Unif' or 'Pois'")
  }
  if ((tolower(ni.type) %in% c("u", "unif", "uniform")) & length(ni) != 2) {
    stop("ni must be a vector of length two if ni.type is uniform ")
  }
  if ((tolower(ni.type) %in% c("m", "max", "maximum")) & length(ni) != 1) {
    stop("ni must be a scalar if ni.type is max ")
  }
  if ((length(ni) == 1) & (ni[1] <= 0)) {
    stop("ni must be strictly positive")
  }
  if ((length(ni) == 2) & (min(ni) == max(ni))) {
    stop("The upper and lower bounds of ni must be different")
  }
  if ((length(ni) == 2) & ((ni[1] < 0) | (ni[2] < 0))) {
    stop("The upper and lower bounds of ni must be strictly positive")
  }
  if ((!tolower(timescale) %in% c("gap", "calendar"))) {
    stop("timescale must be 'gap' or 'calendar'")
  }
  if ((!betaTest.type %in% c("betaRtest", "betaDtest", "joint"))) {
    stop("betaTest.type must be either 'betaRtest', 'betaDtest' or 'joint' ")
  }
  if (ratio < 0) {
    stop("The ratio in favor of the experimental arm must be strictly positive")
  }
  if ((ratio > 4) | (ratio < 0.25)) {
    stop("Such an unbalanced ratio is not recommended")
  }


  P_E <- 1 / (1 / ratio + 1)

  theta <- theta
  eta <- eta
  betaD.HA <- betaD.HA
  betaD.H0 <- betaD.H0
  betaR.HA <- betaR.HA
  betaR.H0 <- betaR.H0

  npar <- 8
  if (test.type == "1-sided") {
    typeIalpha <- typeIerror
  }
  if (test.type == "2-sided") {
    typeIalpha <- typeIerror / 2
  }

  xi_glag <- .glag_xi
  wi_glag <- .glag_wi

  scaleR_W <- medianR.H0 / log(2)^(1 / shapeR.W)
  scaleD_W <- medianD.H0 / log(2)^(1 / shapeD.W)

  h0 <- function(t, y, p) {
    return((t / y)^(p - 1) * (p / y))
  }
  H0 <- function(t, y, p) {
    return((t / y)^p)
  }
  dh0_epsi <- function(t, y, p) {
    return(-p^2 * t^(p - 1) / y^(p + 1))
  }
  dh0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((1 / y) * (t / y)^(p - 1) * (1 + p * log(t / y)))
    }
  }
  dH0_epsi <- function(t, y, p) {
    return(-p * (t^p) / y^(p + 1))
  }
  dH0_nu <- function(t, y, p) {
    if ((t[1] == 0) || (length(t) == 0)) {
      return(0)
    } else {
      return((t / y)^p * log(t / y))
    }
  }

  H0_1 <- function(s, y, p) {
    return(y * s^(1 / p))
  }
  cumsumrev <- function(x) {
    output <- rep(NA, length(x))
    ll <- length(x)
    if (ll > 1) {
      x1 <- x[1:(ll - 1)]
      x2 <- x[2:ll]
      output <- c(x[1], x2 - x1)
    }
    if (ll == 1) {
      output <- x
    }
    return(output)
  }
  G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    return(omega_i^
      {
        Di + deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }

  dbetaR_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, H00_R_Z) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * ((Di * eta + 1) * H00_R_Z / (1 + eta * omega_i * H00_R)^(Di + 1 / eta + 1)))
  }

  dbetaD_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, H00_D_Z) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * H00_D_Z * exp(-omega_i * (1 / theta - 1 + H00_D)) /
      (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }


  dtheta_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    return(omega_i^
      {
        Di + deltaStar + 1 / theta - 1
      } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * omega_i * H00_R)^(Di + 1 / eta) *
      (omega_i - log(omega_i)) / theta^2)
  }

  deta_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D) {
    alpha <- omega_i * H00_R

    deriveta <- exp((Di + 1 / eta) * log(1 + eta * alpha)) * ((-1 / eta^2) * log(1 + eta * alpha) + (Di + 1 / eta) * alpha / (1 + eta * alpha))
    return(-omega_i^{
      Di + deltaStar + 1 / theta - 1
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) / (1 + eta * alpha)^{
      2 * (Di + 1 / eta)
    } *
      deriveta)
  }


  dh0R_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, dH00_R) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta - 1
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * (Di + 1 / eta) * eta * omega_i * dH00_R / (1 + eta * omega_i * H00_R)^(Di + 1 / eta + 1))
  }

  dh0D_G_THETA <- function(omega_i, Di, deltaStar, theta, eta, H00_R, H00_D, dH00_D) {
    return(-omega_i^{
      Di + deltaStar + 1 / theta
    } * exp(-omega_i * (1 / theta - 1 + H00_D)) * dH00_D / (1 + eta * omega_i * H00_R)^(Di + 1 / eta))
  }

  sum_score2_H0 <- sum_score2_HA <- matrix(0, nrow = npar, ncol = npar)
  sum_Di_H0 <- sum_Di_HA <- 0
  Di_H0 <- Di_HA <- 0
  sum_deltaStar_H0 <- sum_deltaStar_HA <- 0


  set.seed(seed)
  for (jt in c(1:samples.mc))
  {
    if (tolower(ni.type) %in% c("u", "unif", "uniform")) {
      if (ni[1] == ni[2]) {
        ni.type <- "max"
        tmp <- ni[1]
        ni <- tmp
        rm(tmp)
      } else {
        ni_ <- round(sample(ni[1]:ni[2], 1))
      }
    }
    if (tolower(ni.type) %in% c("m", "max", "maximum")) {
      ni_ <- round(ni)
    }

    if (tolower(ni.type) %in% c("p", "pois", "poisson")) {
      ni_ <- rpois(n = 1, lambda = ni)
    }

    U_D <- runif(n = 1, 0, 1)

    Xi <- sample(x = c(0, 1), size = 1, prob = c(1 - P_E, P_E))
    Accrual <- runif(n = 1, min = 0, max = Acc.Dur)

    omega_i <- rgamma(1, shape = theta^(-1), scale = theta)
    vi <- rgamma(1, shape = eta^(-1), scale = eta)

    Ti_D_H0 <- H0_1(s = -log(U_D) * exp(-(betaD.H0 * Xi)) / (omega_i), y = scaleD_W, p = shapeD.W)
    Ti_D_HA <- H0_1(s = -log(U_D) * exp(-(betaD.HA * Xi)) / (omega_i), y = scaleD_W, p = shapeD.W)


    if (ni_ > 0) {
      U <- runif(n = ni_, 0, 1)
      Tij_H0 <- H0_1(s = -log(U) * exp(-(betaR.H0 * Xi)) / (omega_i * vi), y = scaleR_W, p = shapeR.W)
      Tij_HA <- H0_1(s = -log(U) * exp(-(betaR.HA * Xi)) / (omega_i * vi), y = scaleR_W, p = shapeR.W)
    } else {
      Tij_H0 <- Ti_D_H0
      Tij_HA <- Ti_D_HA
    }


    if (tolower(FUP.type) %in% c("f", "fixed")) {
      Yi_D_H0 <- pmin(FUP, Ti_D_H0)
      Yi_D_HA <- pmin(FUP, Ti_D_HA)
    }
    if (tolower(FUP.type) %in%  c("uptoend", "end")) {
      Yi_D_H0 <- pmin(Acc.Dur + FUP - Accrual, Ti_D_H0)
      Yi_D_HA <- pmin(Acc.Dur + FUP - Accrual, Ti_D_HA)
    }

    Yij_cal_H0 <- pmin(cumsum(Tij_H0), Yi_D_H0)
    Yij_cal_HA <- pmin(cumsum(Tij_HA), Yi_D_HA)

    deltaij <- as.numeric(Yij_cal_H0 < Yi_D_H0)
    deltaij_ <- as.numeric(Yij_cal_HA < Yi_D_HA)
    Di <- sum(deltaij)
    Di_ <- sum(deltaij_)
    deltaStar <- as.numeric(Ti_D_H0 <= Yi_D_H0)
    deltaStar_ <- as.numeric(Ti_D_HA <= Yi_D_HA)
    DiXi <- sum(deltaij * Xi)
    DiXi_ <- sum(deltaij_ * Xi)

    if (tolower(timescale) == "gap") {
      Yij_H0 <- cumsumrev(Yij_cal_H0)
      Yij_HA <- cumsumrev(Yij_cal_HA)

      R0Yij_H0 <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      R0Yij_HA <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_epsiYij_H0 <- sum(dH0_epsi(Yij_H0, scaleR_W, shapeR.W)) * exp(betaR.H0 * Xi)
      dR0_epsiYij_HA <- sum(dH0_epsi(Yij_HA, scaleR_W, shapeR.W)) * exp(betaR.HA * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      dR0_nuYij_H0 <- sum(dH0_nu(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.H0 * Xi)
      dR0_nuYij_HA <- sum(dH0_nu(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) * exp(betaR.HA * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)

      R0Yij_H0_X <- sum(H0(Yij_H0, scaleR_W, shapeR.W)) * Xi * exp(betaR.H0 * Xi)
      R0Yij_HA_X <- sum(H0(Yij_HA, scaleR_W, shapeR.W)) * Xi * exp(betaR.HA * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    } else {
      Yij_H0 <- Yij_cal_H0
      Yij_HA <- Yij_cal_HA
      n_j_H0 <- length(Yij_cal_H0)
      n_j_HA <- length(Yij_cal_HA)

      if (n_j_H0 > 0) {
        T_max_H0 <- Yij_cal_H0[n_j_H0]
        T_min_H0 <- Yij_cal_H0[1]

        R0Yij_H0 <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_epsiYij_H0 <- (
          dH0_epsi(T_max_H0, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        dR0_nuYij_H0 <- (
          dH0_nu(T_max_H0, scaleR_W, shapeR.W) -
            dH0_nu(T_min_H0, scaleR_W, shapeR.W)
        ) * exp(betaR.H0 * Xi)

        R0Yij_H0_X <- (
          H0(T_max_H0, scaleR_W, shapeR.W) -
            H0(T_min_H0, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.H0 * Xi)
      } else {
        R0Yij_H0 <- 0
        dR0_epsiYij_H0 <- 0
        dR0_nuYij_H0 <- 0
        R0Yij_H0_X <- 0
      }

      H0Yi_H0 <- H0(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_epsiYi_H0 <- dH0_epsi(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      dH0_nuYi_H0 <- dH0_nu(Yi_D_H0, scaleD_W, shapeD.W) * exp(betaD.H0 * Xi)
      H0Yi_H0_X <- H0(Yi_D_H0, scaleD_W, shapeD.W) * Xi * exp(betaD.H0 * Xi)

      if (n_j_HA > 0) {
        T_max_HA <- Yij_cal_HA[n_j_HA]
        T_min_HA <- Yij_cal_HA[1]
        R0Yij_HA <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_epsiYij_HA <- (
          dH0_epsi(T_max_HA, scaleR_W, shapeR.W) -
            dH0_epsi(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        dR0_nuYij_HA <- (
          dH0_nu(T_max_HA, scaleR_W, shapeR.W) -
            dH0_nu(T_min_HA, scaleR_W, shapeR.W)
        ) * exp(betaR.HA * Xi)

        R0Yij_HA_X <- (
          H0(T_max_HA, scaleR_W, shapeR.W) -
            H0(T_min_HA, scaleR_W, shapeR.W)
        ) * Xi * exp(betaR.HA * Xi)
      } else {
        R0Yij_HA <- 0
        dR0_epsiYij_HA <- 0
        dR0_nuYij_HA <- 0
        R0Yij_HA_X <- 0
      }

      H0Yi_HA <- H0(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_epsiYi_HA <- dH0_epsi(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      dH0_nuYi_HA <- dH0_nu(Yi_D_HA, scaleD_W, shapeD.W) * exp(betaD.HA * Xi)
      H0Yi_HA_X <- H0(Yi_D_HA, scaleD_W, shapeD.W) * Xi * exp(betaD.HA * Xi)
    }

    G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    dbetaR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_R_Z = R0Yij_H0_X)
    }))
    dbetaR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaR_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_R_Z = R0Yij_HA_X)
    }))
    score_betaR_H0 <- DiXi + dbetaR_G_THETA_val_H0 / G_THETA_val_H0
    score_betaR_HA <- DiXi_ + dbetaR_G_THETA_val_HA / G_THETA_val_HA

    dbetaD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, H00_D_Z = H0Yi_H0_X)
    }))
    dbetaD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dbetaD_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, H00_D_Z = H0Yi_HA_X)
    }))
    score_betaD_H0 <- deltaStar * Xi + dbetaD_G_THETA_val_H0 / G_THETA_val_H0
    score_betaD_HA <- deltaStar_ * Xi + dbetaD_G_THETA_val_HA / G_THETA_val_HA

    dtheta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    dtheta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dtheta_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))

    score_theta_H0 <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_H0 / G_THETA_val_H0
    score_theta_HA <- digamma(1 / theta) / theta^2 + (log(theta) - 1) / theta^2 + dtheta_G_THETA_val_HA / G_THETA_val_HA

    deta_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0)
    }))
    deta_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      deta_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA)
    }))
    l <- l_ <- 0
    if (Di > 0) {
      l <- 1:Di
    }
    if (Di_ > 0) {
      l_ <- 1:Di_
    }

    score_eta_H0 <- Di / eta - sum((1 / eta^2) / (Di + 1 / eta - l)) + deta_G_THETA_val_H0 / G_THETA_val_H0
    score_eta_HA <- Di_ / eta - sum((1 / eta^2) / (Di_ + 1 / eta - l_)) + deta_G_THETA_val_HA / G_THETA_val_HA

    depsiR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_epsiYij_H0)
    }))
    depsiR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_epsiYij_HA)
    }))

    score_epsiR_H0 <- -Di * shapeR.W / scaleR_W + depsiR_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiR_HA <- -Di_ * shapeR.W / scaleR_W + depsiR_G_THETA_val_HA / G_THETA_val_HA


    dnuR_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_R = dR0_nuYij_H0)
    }))
    dnuR_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0R_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_R = dR0_nuYij_HA)
    }))

    score_nuR_H0 <- sum(deltaij * dh0_nu(Yij_H0, scaleR_W, shapeR.W) / h0(Yij_H0, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_H0 / G_THETA_val_H0
    score_nuR_HA <- sum(deltaij_ * dh0_nu(Yij_HA, scaleR_W, shapeR.W) / h0(Yij_HA, scaleR_W, shapeR.W), na.rm = TRUE) + dnuR_G_THETA_val_HA / G_THETA_val_HA

    depsiD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_epsiYi_H0)
    }))
    depsiD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_epsiYi_HA)
    }))

    score_epsiD_H0 <- -deltaStar * shapeD.W / scaleD_W + depsiD_G_THETA_val_H0 / G_THETA_val_H0
    score_epsiD_HA <- -deltaStar_ * shapeD.W / scaleD_W + depsiD_G_THETA_val_HA / G_THETA_val_HA

    dnuD_G_THETA_val_H0 <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di, deltaStar = deltaStar, theta = theta, eta = eta, H00_R = R0Yij_H0, H00_D = H0Yi_H0, dH00_D = dH0_nuYi_H0)
    }))
    dnuD_G_THETA_val_HA <- sum(wi_glag * sapply(xi_glag, function(x) {
      dh0D_G_THETA(omega_i = x, Di = Di_, deltaStar = deltaStar_, theta = theta, eta = eta, H00_R = R0Yij_HA, H00_D = H0Yi_HA, dH00_D = dH0_nuYi_HA)
    }))

    score_nuD_H0 <- deltaStar * dh0_nu(Yi_D_H0, scaleD_W, shapeD.W) / h0(Yi_D_H0, scaleD_W, shapeD.W) + dnuD_G_THETA_val_H0 / G_THETA_val_H0
    score_nuD_HA <- deltaStar_ * dh0_nu(Yi_D_HA, scaleD_W, shapeD.W) / h0(Yi_D_HA, scaleD_W, shapeD.W) + dnuD_G_THETA_val_HA / G_THETA_val_HA

    if (npar == 8) {
      score_vector_H0 <- c(score_betaR_H0, score_betaD_H0, score_theta_H0, score_eta_H0, score_epsiR_H0, score_nuR_H0, score_epsiD_H0, score_nuD_H0)
      score_vector_HA <- c(score_betaR_HA, score_betaD_HA, score_theta_HA, score_eta_HA, score_epsiR_HA, score_nuR_HA, score_epsiD_HA, score_nuD_HA)
    }

    U_score2_H0 <- (score_vector_H0) %*% t(score_vector_H0)
    U_score2_HA <- (score_vector_HA) %*% t(score_vector_HA)


    sum_score2_H0 <- sum_score2_H0 + U_score2_H0
    sum_score2_HA <- sum_score2_HA + U_score2_HA

    if (ni_ > 0) {
      sum_Di_H0 <- sum_Di_H0 + Di / ni_
      Di_H0 <- Di_H0 + Di
      sum_Di_HA <- sum_Di_HA + Di_ / ni_
      Di_HA <- Di_HA + Di_
    }
    sum_deltaStar_H0 <- sum_deltaStar_H0 + deltaStar
    sum_deltaStar_HA <- sum_deltaStar_HA + deltaStar_
  }
  E_score2_H0 <- sum_score2_H0 / samples.mc
  E_score2_HA <- sum_score2_HA / samples.mc
  E_Di_H0 <- sum_Di_H0 / samples.mc
  E_Di_HA <- sum_Di_HA / samples.mc
  E_deltaStar_H0 <- sum_deltaStar_H0 / samples.mc
  E_deltaStar_HA <- sum_deltaStar_HA / samples.mc


  if (is.square.matrix(E_score2_H0) & is.positive.definite(E_score2_H0)) {
    E1_H_H0 <- chol2inv(chol(E_score2_H0))
  } else {
    E1_H_H0 <- solve(E_score2_H0)
  }
  if (is.square.matrix(E_score2_HA) & is.positive.definite(E_score2_HA)) {
    E1_H_HA <- chol2inv(chol(E_score2_HA))
  } else {
    E1_H_HA <- solve(E_score2_HA)
  }


  C_trt <- switch(betaTest.type,
    "betaRtest" = matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "betaDtest" = matrix(data = c(0, 1, 0, 0, 0, 0, 0, 0), nrow = 1, ncol = npar),
    "joint" = matrix(data = c(
      1, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 2, ncol = npar)
  )

  cat("\n")


  if (statistic == "Wald") {
    rr <- nrow(C_trt)
    ncp <- .chi2ncp(typeIalpha * 2, power, df = rr)
    ncp0 <- 0
    CIC <- t(C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)))

    Npts <- ceiling(ncp * solve(CIC))

    ncpA <- (Npts) * t(C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0))) %*% solve(C_trt %*% E1_H_HA %*% t(C_trt)) %*% (C_trt %*% (c(betaR.HA, betaD.HA, 0, 0, 0, 0, 0, 0)))
    estimated.power <- 1 - pchisq(qchisq(1 - typeIalpha * 2, df = rr, ncp = ncp0), df = rr, ncp = ncpA)
  }

  res <- list(
    Npts = Npts,
    estimated.power = estimated.power,
    events.rec = c(Npts * Di_H0 / samples.mc, Npts * Di_HA / samples.mc),
    events.D = c(Npts * sum_deltaStar_H0 / samples.mc, Npts * sum_deltaStar_HA / samples.mc)
    # censoring.rec = c((1 - E_Di_H0) * 100, (1 - E_Di_HA) * 100),
    # censoring.D = c((1 - E_deltaStar_H0) * 100, (1 - E_deltaStar_HA) * 100)
  )

  class(res) <- "frailtyDesign"

  return(res)
}
