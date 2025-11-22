#' Plot Method for a Weibull competing risks model with optional shared frailty between transitions.
#'
#' @name plot.frailtyCmprsk
#' @aliases  plot.frailtyCmprsk
#' @title Plot Method for a Weibull competing risks model with optional shared frailty between transitions.
#'
#' @description
#' Plots estimated baseline survival and hazard functions from an object of
#' class 'frailtyCmprsk'. Confidence bands are allowed.
#' 
#' 


#' @usage
#' \method{plot}{frailtyCmprsk}(x, type.plot="Baseline hazard", events, conf.bands=TRUE, 
#' pos.legend="topright", cex.legend=0.7, lwd=c(1,1,1), color=2, median=TRUE, 
#' Xlab = "Time", Ylab = "Baseline hazard function",...)
#'
#' @param x A Weibull competing risks model, i.e. an \code{frailtyCmprsk} class object
#'   (output from calling \code{frailtyCmprsk} function).
#' @param type.plot a character string specifying the type of curve. Possible
#'   value are "Baseline hazard", or "Baseline survival". The default is "Baseline hazard".
#'   Only the first letters are required, e.g "Haz", "Su".
#' @param events Integer vector specifying which competing events to display in the plots. 
#'   Use 1 for the first event, 2 for the second, and so on. If not specified, plots for all events are shown.
#' @param conf.bands Logical value. Determines whether confidence bands will be
#'   plotted.  The default is TRUE.
#' @param pos.legend The location of the legend can be specified by setting
#'   this argument to a single keyword from the list '"bottomright"', '"bottom"',
#'   '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#'   '"center"'. The default is '"topright"'
#' @param cex.legend character expansion factor *relative* to current
#'   'par("cex")'. Default is 0.5
#' @param lwd A vector of length 3 of positive values to specify the line width of the plots and their
#'   confidence bands. If no confidence bands are plotted, the last two elements of the vector are ignored.
#'   Default is (1,1,1).
#' @param color color of the curve (integer)
#' @param median Logical value. Determines whether median survival time will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'
#' @param Ylab Label of y-axis. Default is '"Baseline hazard function"'
#' @param ... other unused arguments

#'
#' @return
#' Print a plot of a Weibull competing risks model with optional shared
#'  frailty between transitions.
#'
#' @seealso
#'   \code{\link{frailtyCmprsk}}
#'
#' @keywords file
#'
#' @examples
#' \donttest{
#'
#'     ###--- Simple Weibull competing risks model ---###
#'
#'     data(CPRSKbmtcrr)
#'
#'     ModCmprsk_NoCov_Factor <- frailtyCmprsk(
#'       formulas = list(
#'         Surv(observed_time, Status, type = "mstate") ~ D,
#'         ~ 1
#'       ),
#'       data = CPRSKbmtcrr,
#'       print.info = FALSE,
#'       maxit = 100
#'     )
#'
#'     plot(ModCmprsk_NoCov_Factor)
#'
#'     #-- No confidence bands
#'     plot(ModCmprsk_NoCov_Factor, conf.bands = FALSE)
#'
#' }
#'
#' @import graphics 






 
#' @export
"plot.frailtyCmprsk" <- function (x, type.plot="Baseline hazard", events, conf.bands=TRUE, pos.legend="topright", 
                                cex.legend=0.7, lwd=c(1,1,1), color=2, median=TRUE, Xlab = "Time",
                                Ylab = "Baseline hazard function",...)
{
  
  
  
  
  # Function to calculate confidence bands for a set of times
  
  hazard_and_confidence_bands_frailtycmprsk <- function(t_vec,scale,shape,varcov) {
    
    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    
    ##DELTA METHOD TO OBTAIN THE VARCOV MATRIX FOR LOG OF SHAPE AND SCALE
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)
    
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline hazard estimate","Lower", "Upper")    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      haz_pert <- numeric(2000)
      for (k in 1:2000) {
        haz_pert[k] <- (Sample_scale_shape[k,2]/ (Sample_scale_shape[k,1]^Sample_scale_shape[k,2])) * (t^(Sample_scale_shape[k,2] - 1))
      }
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      result[t_idx, 3] <- quantile(haz_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(haz_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  survival_and_confidence_bands_frailtycmprsk <- function(t_vec,scale,shape,varcov) {
    
    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)  
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline survival estimate","Lower", "Upper")
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      surv_pert <- numeric(2000)
      for (k in 1:2000) {
        surv_pert[k] <- exp(-(t/Sample_scale_shape[k,1])^Sample_scale_shape[k, 2])  }
      
      result[t_idx,1] <- t
      result[t_idx,2] <- exp(-(t/scale)^shape)
      
      result[t_idx, 3] <- quantile(surv_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(surv_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  base_hazard_frailtycmprsk <- function(t_vec,scale,shape) {
    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline hazard estimate")
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      
    }
    
    return(result)
  }
  
  
  
  
  
  su_frailtycmprsk <- function(t_vec,scale,shape) {
    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline survival estimate")
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      result[t_idx,1] <- t
      result[t_idx,2] <- exp(-(t/scale)^shape)
      
      
    }
    
    return(result)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(!missing(lwd))
    if(!all(sapply(lwd, function(x) is.numeric(x) && x>=0  ))) {
      stop("lwd should be a vector of positive real numbers.")
    }
  
  if(length(lwd)!=3){
    stop("'lwd should be a vector of size 3'.")
  }
  
  
  
  if (!inherits(x, "frailtyCmprsk")) {
    stop("'x' must be an 'frailtyCmprsk' model object")
  }
  
  
  vcov <- x$vcov
  partialH <- x$partialH
  data <- x$data
  
  y=data$y
  delta=data$delta
  n_competing_events_from_formulas <- length(x$formulas)
  delta_indicators <- lapply(1:n_competing_events_from_formulas, function(k) ifelse(delta == k, 1, 0))
  
  
  
  x0 <- vector("list", n_competing_events_from_formulas)
  for (k in 1:n_competing_events_from_formulas) {
    
    max_time_for_event_k <- max(y[delta_indicators[[k]] == 1])
    x0[[k]] <- round(seq(0, max_time_for_event_k, length.out = 99), 3)
  }
  
  plot.type <- charmatch(type.plot, c("Baseline hazard", "Haz", "Baseline survival", "Su"), nomatch = 0)
  
  
  
 
  
  
  if (plot.type == 0) {
    stop("estimator must be 'Baseline hazard' (or 'Haz'), 'Baseline survival' (or 'Su')")
  }
  
  
  if(plot.type==3 || plot.type==4){
    if((median)){
      median=TRUE
    }
    
  }
  if(!(plot.type==3 || plot.type==4)){
    
    median=FALSE
    
  }
  
  if(!missing(conf.bands)){
    if (!is.logical(conf.bands) || length(conf.bands) != 1) {
      stop("'conf.bands' should be a logical value (TRUE or FALSE).")
    }
  }
  
  
  if(!missing(median)){
    if (!is.logical(median) || length(median) != 1) {
      stop("'median' should be a logical value (TRUE or FALSE).")
    }
  }
  
  
  if(!missing(cex.legend)){
    if (!is.numeric(cex.legend) || length(cex.legend) != 1 || cex.legend <= 0) {
      stop("'cex.legend' should be a positive numeric value.")
    }
  }
  
 
  
  if (!missing(color) && (!is.numeric(color) || length(color) != 1 || color <= 0 || color != round(color))) {
    stop("'color' should be a positive integer.")
  }
  
  
  
  valid_pos_legend <- c("topright","topleft","bottomright","bottomleft","top","bottom","left","right","center")
  if(!missing(pos.legend)){
    if (!(pos.legend %in% valid_pos_legend)) {
      stop("'pos.legend' should be one of: 'topright','topleft','bottomright','bottomleft','top','bottom','left','right','center'.")
    }
    
  }
  
  
  
  valid_event_indices <- 1:length(x$formulas) 
  if (!missing(events)) {
    events <- unique(events)
    if (!is.numeric(events) || !all(events == round(events))) {
      stop("'events' must be a vector of integers.")
    }
    
    if (!all(events %in% valid_event_indices)) {
      valid_events_str <- paste(valid_event_indices, collapse = ", ")
      stop(paste0("All values in 'events' must be integers corresponding to ",
                  "competing events (1 to ", length(x$formulas), "). ",
                  "Valid event indices are: ", valid_events_str, "."))
    }
    
    
  }
  
  
 
  
  if(!missing(Xlab)){
    if (!is.character(Xlab) || length(Xlab) != 1) {
      stop("'Xlab' should be a single character string.")
    }
  }
  
  if(!missing(Ylab)){
    if (!is.character(Ylab) || length(Ylab) != 1) {
      stop("'Ylab' should be a single character string.")
    }
  }
  
  
  
  model <- x
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar))
  
  
  
  ### HAZARD
  if(!missing(events)){
    if(plot.type==1 || plot.type==2){
      if(conf.bands==TRUE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          if(missing(Ylab)){
            Ylab <- paste0("Baseline hazard function")
          }
          no_conf <- NULL
          if((2*k-1) %in% partialH | (2*k) %in% partialH){
            no_conf <- 1
            message(paste0("Some parameters of the baseline hazard for event ",k," were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case."))
            lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!((2*k-1) %in% partialH | (2*k) %in% partialH)) && dim(model$lam0[[k]])[2]==2 ){
            no_conf <- 1
            message(paste0("Numerical problem encountered in computing the confidence bands of the baseline hazard for event ",k))
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            lam0 <- hazard_and_confidence_bands_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k],vcov[(2*k-1):(2*k), (2*k-1):(2*k)])
            matplot(x0[[k]], lam0[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
      }
      
      
      
      
      if(conf.bands==FALSE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          
          if(missing(Ylab)){
            Ylab <- paste0("Baseline hazard function")
          }
          
          lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
          matplot(x0[[k]], lam0[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = 1, cex = cex.legend)
          
          
          
        }
        
        
        
      }
      
      
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(missing(events)){
    events <- c(1:n_competing_events_from_formulas)
    if(plot.type==1 || plot.type==2){
      if(conf.bands==TRUE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          if(missing(Ylab)){
            Ylab <- paste0("Baseline hazard function")
          }
          no_conf <- NULL
          if((2*k-1) %in% partialH | (2*k) %in% partialH){
            no_conf <- 1
            message(paste0("Some parameters of the baseline hazard for event ",k," were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case."))
            lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!((2*k-1) %in% partialH | (2*k) %in% partialH)) && dim(model$lam0[[k]])[2]==2 ){
            no_conf <- 1
            message(paste0("Numerical problem encountered in computing the confidence bands of the baseline hazard for event ",k))
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            lam0 <- hazard_and_confidence_bands_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k],vcov[(2*k-1):(2*k), (2*k-1):(2*k)])
            matplot(x0[[k]], lam0[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
      }
      
      
      
      
      if(conf.bands==FALSE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          
          if(missing(Ylab)){
            Ylab <- paste0("Baseline hazard function")
          }
          
          lam0 <- base_hazard_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
          matplot(x0[[k]], lam0[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = 1, cex = cex.legend)
          
          
          
        }
        
        
        
      }
      
      
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### SURVIVAL
  
  
  
  if(!missing(events)){
    if(plot.type==3 || plot.type==4){
      if(conf.bands==TRUE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          if(missing(Ylab)){
            Ylab <- paste0("Baseline survival function")
          }
          no_conf <- NULL
          if((2*k-1) %in% partialH | (2*k) %in% partialH){
            no_conf <- 1
            message(paste0("Some parameters of the baseline survival for event ",k," were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case."))
            surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], surv0[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
          
          if((!((2*k-1) %in% partialH | (2*k) %in% partialH)) && dim(model$surv0[[k]])[2]==2 ){
            no_conf <- 1
            message(paste0("Numerical problem encountered in computing the confidence bands of the baseline survival for event ",k))
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
          
          
          if(is.null(no_conf)){
            surv0 <- survival_and_confidence_bands_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k],vcov[(2*k-1):(2*k), (2*k-1):(2*k)])
            matplot(x0[[k]], surv0[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
        }
      }
      
      
      
      
      if(conf.bands==FALSE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          
          if(missing(Ylab)){
            Ylab <- paste0("Baseline survival function")
          }
          
          surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
          matplot(x0[[k]], surv0[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
          if (median == TRUE) { 
            
            abline(h = 0.5, col = "black", lty = 2, lwd = 1)
            
          }
          
          
          legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
          legend_col <- if (is.null(no_conf)) rep(color, 2) else color
          legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
          
          if (median == TRUE ) {
            legend_text <- c(legend_text, "Median survival")
            legend_col <- c(legend_col, "black")
            legend_lty <- c(legend_lty, 2)
          }
          
          legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
          
          
          
        }
        
        
        
      }
      
      
      
    }
  }
  
  
  
  
  
  
  
  
  if(missing(events)){
    events <- c(1:n_competing_events_from_formulas)
    if(plot.type==3 || plot.type==4){
      if(conf.bands==TRUE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          if(missing(Ylab)){
            Ylab <- paste0("Baseline survival function")
          }
          no_conf <- NULL
          if((2*k-1) %in% partialH | (2*k) %in% partialH){
            no_conf <- 1
            message(paste0("Some parameters of the baseline survival for event ",k," were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case."))
            surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], surv0[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
          
          if((!((2*k-1) %in% partialH | (2*k) %in% partialH)) && dim(model$surv0[[k]])[2]==2 ){
            no_conf <- 1
            message(paste0("Numerical problem encountered in computing the confidence bands of the baseline survival for event ",k))
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
            matplot(x0[[k]], lam0[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
          
          
          if(is.null(no_conf)){
            surv0 <- survival_and_confidence_bands_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k],vcov[(2*k-1):(2*k), (2*k-1):(2*k)])
            matplot(x0[[k]], surv0[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
            if (median == TRUE) { 
              
              abline(h = 0.5, col = "black", lty = 2, lwd = 1)
              
            }
            
            
            legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
            legend_col <- if (is.null(no_conf)) rep(color, 2) else color
            legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
            
            if (median == TRUE ) {
              legend_text <- c(legend_text, "Median survival")
              legend_col <- c(legend_col, "black")
              legend_lty <- c(legend_lty, 2)
            }
            
            legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
            
            
          }
        }
      }
      
      
      
      
      if(conf.bands==FALSE){
        par(mfrow=c(1,length(events)))
        for(k in events){
          
          if(missing(Ylab)){
            Ylab <- paste0("Baseline survival function")
          }
          
          surv0 <- su_frailtycmprsk(x0[[k]],x$scale.weib[k],x$shape.weib[k])
          matplot(x0[[k]], surv0[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main=paste0("Event ", k),lwd=lwd)
          if (median == TRUE) { 
            
            abline(h = 0.5, col = "black", lty = 2, lwd = 1)
            
          }
          
          
          legend_text <- if (is.null(no_conf)) c("Baseline survival", "Confidence bands") else c("Baseline survival")
          legend_col <- if (is.null(no_conf)) rep(color, 2) else color
          legend_lty <- if (is.null(no_conf)) c(1, 2) else 1
          
          if (median == TRUE ) {
            legend_text <- c(legend_text, "Median survival")
            legend_col <- c(legend_col, "black")
            legend_lty <- c(legend_lty, 2)
          }
          
          legend(pos.legend, legend = legend_text, col = legend_col, lty = legend_lty, cex = cex.legend)
          
          
          
        }
        
        
        
      }
      
      
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  return(invisible())
}






