#' Generating from a joint competing Joint frailty model with a recurrent event and two
#' terminal events.
#' 
#' @description Generates data under a joint frailty model for a single recurrent event
#' and two terminal events in a calendar-time format. Only a single covariate representing the treatment
#' is allowed. Event times are generated under Weibull distributions. 
#' @usage 
#' simulatejointRecCompet(n, censoring = 28, maxrecurrent = 50, 
#'            par0 = c(shapeR = 1.5, scaleR = 10, 
#'            shapeM = 1.75, scaleM = 16, 
#'            shapeD = 1.75, scaleD = 16, 
#'            sigma = 0.5, alphaM = 1, alphaD = 0, 
#'            betaR = -0.5, betaM = -0.5, betaD = 0))
#' @param n Number of subjects. Default is 1500. 
#' @param censoring A number indicating a fixed right censoring time for all subjects 
#' (as an administrative censoring). If \code{NULL}, no censoring is applied. 
#' Default is 28. 
#' @param maxrecurrent Maximum number of recurrent events per subject.
#' If \code{NULL}, no upper bound is set for the number of of recurrent events
#' Default is 50. 
#' @param par0 A vector of arguments controlling the parameters of the generating model. 
#'  \itemize{
#'   \item{\code{shapeR}: \code{shape} parameter of the Weibull distribution for the recurrent event}
#'   \item{\code{scaleR}: \code{scale} parameter of the Weibull distribution for the recurrent event} 
#'   \item{\code{shapeT1}: \code{shape} parameter of the Weibull distribution for the first terminal event}
#'   \item{\code{scaleT1}: \code{scale} parameter of the Weibull distribution for the first terminal event}
#'   \item{\code{shapeT2}: \code{shape} parameter of the Weibull distribution for the second terminal event}
#'   \item{\code{scaleT2}: \code{scale} parameter of the Weibull distribution for the second terminal event}
#'   \item{\code{sigma}: Standard deviation of the frailty}
#'   \item{\code{alphaT1}: Power parameter (link) of the frailty for the first terminal event}
#'   \item{\code{alphaT2}: Power parameter (link) of the frailty for the second terminal event}
#'   \item{\code{betaR}: Association parameter for the treatment effect on the recurrent event}
#'   \item{\code{betaT1}: Association parameter for the treatment effect on the first terminal event}
#'   \item{\code{betaT2}: Association parameter for the treatment effect on the second terminal event}
#' }
#' @return  Returns a \code{data.frame} object with the following columns:
#' \itemize{
#'  \item{\code{id} Id number for each subject}
#'  \item{\code{treatment} Binary treatment indicator} 
#'  \item{\code{tstart} Start time of the observation period}
#'  \item{\code{tstop} Stop time of the observation period}
#'  \item{\code{recurrent} Censoring indicator for the recurrent event}
#'  \item{\code{terminal1} Censoring indicator for the first terminal event}
#'  \item{\code{terminal2} Censoring indicator for the second terminal event}
#' }
#' @seealso \link{jointRecCompet}
#' @importFrom stats rweibull pweibull dweibull qweibull rbinom rnorm
#' @aliases simulatejointRecCompet 
#' @export
"simulatejointRecCompet"<-function (n = 1500, censoring = 28, maxrecurrent = 50, par0 = c(shapeR = 1.5, scaleR = 10, 
    shapeM = 1.75, scaleM = 16, shapeD = 1.75, scaleD = 16, sigma = 0.5, 
    alphaM = 1, alphaD = 0, betaR = -0.5, betaM = -0.5, betaD = 0)) 
{   

    rweibRH <- function(n, shape ,scale , rh){
                rweibull(n, shape = shape, scale = scale * rh^(-1/shape))
            }

    pweibRH <- function(q, shape, scale, rh){
                pweibull(q, shape = shape, scale = scale * rh^(-1/shape))
            }

    dweibRH <- function(x, shape, scale, rh){
                dweibull(x, shape = shape, scale = scale * rh^(-1/shape))
            }

    qweibRH <- function(p, shape, scale, rh){
                qweibull(p, shape = shape, scale = scale * rh^(-1/shape))
            }
    K = censoring
    trt = rbinom(n, 1, 0.5)
    w = rnorm(n, 0,par0["sigma"])
    T1 = rweibRH(n, shape = par0["shapeM"], 
                scale = par0["scaleM"], 
                rh = exp(w * par0["alphaM"] + 
                            par0["betaM"] * trt))
    T2 = rweibRH(n, shape = par0["shapeD"], 
                scale = par0["scaleD"], 
                rh = exp(w * par0["alphaD"] + 
                            par0["betaD"] * trt))
    y = pmin(T1, T2, K)
    terminal1 = as.numeric(T1 < T2 & T1 < K)
    terminal2 = as.numeric(T2 < T1 & T2 < K)
    t = lapply(seq_along(trt),function(i){
        cumsum(rweibRH(maxrecurrent, shape = par0["shapeR"], 
                    scale = par0["scaleR"], rh = exp(w[i] + par0["betaR"] * 
                                                        trt[i])))
    })
    
    data= do.call(rbind,lapply(1:n,function(id){
        datatemp=data.frame(do.call(rbind,lapply(1:50,function(i){
        trt=trt[id]
        tstart=c(0,t[[id]][1:(length(t)-1)])
        c(id=id,trt=trt,tstart=tstart[i],t=t[[id]][i],term1=terminal1[id],term2=terminal2[id],y=y[id])
        })))
        datatemp=datatemp[datatemp$tstart<datatemp$y,]
        datatemp$terminal1=datatemp$term1*(datatemp$t>datatemp$y)
        datatemp$terminal2=datatemp$term2*(datatemp$t>datatemp$y)
        datatemp$event = 1*(datatemp$t<datatemp$y)
        datatemp$t = sapply(seq_along(datatemp$t),function(i){
        ifelse(datatemp$t[i]>datatemp$y[i],datatemp$y[i],datatemp$t[i])})
        datatemp$y <- datatemp$term1 <- datatemp$term2 <- NULL
        datatemp
    }))
    data=as.data.frame(data)
    data=data[,c(1,2,3,4,7,5,6)]
    colnames(data) = c("id","treatment","tstart","tstop","event","terminal1","terminal2")
    data
}