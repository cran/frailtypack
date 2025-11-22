#' @name plot.frailtyIllnessDeath

#' @aliases plot.frailtyIllnessDeath
#' @title Plot Method for a Weibull Illness-Death model with optional shared frailty between transitions.
#' @usage
#'   \method{plot}{frailtyIllnessDeath}(x, type.plot = "Baseline hazard", transition,
#'    conf.bands=TRUE, pos.legend = "topright", cex.legend=0.7, lwd=c(1,1,1), 
#'    color=2, median=TRUE, Xlab = "Time", Ylab= "Baseline hazard function",...)
#' @param x A Weibull Illness-Death model, i.e. an \code{IllnessDeath} class object
#'     (output from calling \code{IllnessDeath} function).
#' @param type.plot a character string specifying the type of curve. Possible
#'     value are "Baseline hazard", or "Baseline survival". The default is "Baseline hazard".
#'     Only the first letters are required, e.g "Haz", "Su".
#' @param transition Argument to specify if only the plot of one of the three transitions is wanted. If not specified, the plots
#'   for all transitions are displayed. Possible values are "01", "02" or "12".
#' @param conf.bands Logical value. Determines whether confidence bands will be
#'     plotted.  The default is to do so.
#' @param pos.legend The location of the legend can be specified by setting
#'     this argument to a single keyword from the list '"bottomright"', '"bottom"',
#'     '"bottomleft"', '"left"', '"topleft"', '"top"', '"topright"', '"right"' and
#'     '"center"'. The default is '"topright"'.
#' @param cex.legend character expansion factor *relative* to current
#'     'par("cex")'. Default is 0.7.
#' @param lwd A vector of length 3 of positive values to specify the line width of the plots and their
#'   confidence bands. If no confidence bands are plotted, the last two elements of the vector are ignored.
#'   Default is (1,1,1).     
#' @param color color of the curve (integer).
#' @param median Logical value. Determines whether median survival time will be plotted. Default is TRUE.
#' @param Xlab Label of x-axis. Default is '"Time"'.
#' @param Ylab Label of y-axis. Default is '"Baseline hazard function"'.
#' @param ... other unused arguments
#' 
#' @return
#'   Print a plot of a Weibull Illness-Death model with optional shared frailty between transitions.
#' @description
#'   Plots estimated baseline survival and hazard functions from an object of
#'   class 'frailtyIllnessDeath'. Confidence bands are allowed.
#' @examples
#'
#'
#'   \donttest{
#'
#'     ###--- Semi-Markovian Weibull Illness-Death model with left truncation ---###
#'
#'     data(Paq810)
#'
#'     ModIllnessDeath_SemiMarkov_LeftTrunc <- frailtyIllnessDeath(
#'       formula = Surv(e, r, dementia) ~ gender + certif,
#'       formula.terminalEvent = Surv(t, death) ~ gender + certif,
#'       data = Paq810,
#'       print.info = FALSE,
#'       maxit = 100
#'     )
#'
#'     plot(ModIllnessDeath_SemiMarkov_LeftTrunc)
#'
#'     #-- No confidence bands
#'     plot(ModIllnessDeath_SemiMarkov_LeftTrunc, conf.bands = FALSE)
#'
#'   }
#'
#'
#' @seealso
#'   \code{\link{frailtyIllnessDeath}}
#' @keywords file

#' @import graphics 








#' @export
"plot.frailtyIllnessDeath" <- function (x, type.plot="Baseline hazard", transition, conf.bands=TRUE, pos.legend="topright",
                               cex.legend=0.7,lwd=c(1,1,1), color=2, median=TRUE, 
                               Xlab = "Time", Ylab = "Baseline hazard function",...
                               )
{
  
  
  
  
  # Function to calculate confidence bands for a set of times
  hazard_and_confidence_bands_ID <- function(t_vec,scale,shape,varcov) {
    # Extract the scale and shape Weibull parameters
    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    
    ##DELTA METHOD TO OBTAIN THE VARCOV MATRIX FOR LOG OF SHAPE AND SCALE
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)
    
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    # Initialize a matrix to store hazard quantiles
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline hazard estimate","Lower", "Upper")    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      haz_pert <- numeric(2000)
      for (k in 1:2000) {
        haz_pert[k] <- (Sample_scale_shape[k,2]/ (Sample_scale_shape[k,1]^Sample_scale_shape[k,2])) * (t^(Sample_scale_shape[k,2] - 1))
      }
      # Compute the quantiles for the current time
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      result[t_idx, 3] <- quantile(haz_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(haz_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  # Function to calculate confidence bands for a set of times
  survival_and_confidence_bands_ID <- function(t_vec,scale,shape,varcov) {
    # Extract the scale and shape Weibull parameters
    scale <- scale
    shape <- shape
    log_scale <- log(scale)
    log_shape <- log(shape)
    b <- c(log_shape, log_scale)
    ##DELTA METHOD TO OBTAIN THE VARCOV MATRIX FOR SQRT OF SHAPE AND SCALE
    varcov <- diag(c(1/scale,1/shape),2,2) %*% varcov %*% diag(c(1/scale,1/shape),2,2)  # Perform Cholesky decomposition of the covariance matrix
    
    
    # Initialize a matrix to store hazard quantiles
    result <- matrix(0, nrow = length(t_vec), ncol = 4)
    colnames(result) <- c("Time","Baseline survival estimate","Lower", "Upper")
    
    sample_logscale_logshape <-  mvrnorm(n = 2000, mu = c(log_scale, log_shape), Sigma = varcov)
    Sample_scale_shape <- exp(sample_logscale_logshape)
    
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      surv_pert <- numeric(2000)
      for (k in 1:2000) {
        surv_pert[k] <- exp(-(t/Sample_scale_shape[k,1])^Sample_scale_shape[k, 2])  }
      # Compute the quantiles for the current time
      # Compute the quantiles for the current time
      result[t_idx,1] <- t
      result[t_idx,2] <- exp(-(t/scale)^shape)
      
      result[t_idx, 3] <- quantile(surv_pert, prob = 0.025,na.rm=TRUE)
      result[t_idx, 4] <- quantile(surv_pert, prob = 0.975,na.rm=TRUE)
    }
    
    return(result)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  # Function to calculate hazars for a set of times
  base_hazard_ID <- function(t_vec,scale,shape) {
    # Extract the scale and shape Weibull parameters
    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline hazard estimate")
    
    # Loop over each time in t_vec
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      # Compute the quantiles for the current time
      result[t_idx,1] <- t
      result[t_idx,2] <- (shape/ (scale^shape)) * (t^(shape - 1))
      
      
    }
    
    return(result)
  }
  
  
  
  
  
  # Function to calculate hazars for a set of times
  su_ID <- function(t_vec,scale,shape) {
    # Extract the scale and shape Weibull parameters
    scale <- scale
    shape <- shape
    
    
    
    result <- matrix(0, nrow = length(t_vec), ncol = 2)
    colnames(result) <- c("Time","Baseline survival estimate")
    
    # Loop over each time in t_vec
    for (t_idx in seq_along(t_vec)) {
      t <- t_vec[t_idx]
      
      
      # Compute the quantiles for the current time
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
  
  
  
  if (!inherits(x, "frailtyIllnessDeath")) {
    stop("'x' must be an 'frailtyIllnessDeath' model object")
  }
  
  
  vcov <- x$vcov
  partialH <- x$partialH
  data <- x$data
  l=data$l
  y1=data$y1
  y2=data$y2
  delta1=data$delta1
  delta2=data$delta2
  
  
  x01=round(seq(min(l),max(y1),length=99),3)
  x02=round(seq(min(l),max(y2[which(delta1==0)]),length=99),3)
  
  if(x$model =="Semi-Markov"){
    
  x12=round(seq(0,max(y2[which(delta1==1)]-y1[which(delta1==1)]),length=99),3)
    }else{
      x12=round(seq(min(y1[which(delta1==1)]),max(y2[which(delta1==1)]),length=99),3)
    }

  
  plot.type <- charmatch(type.plot, c("Baseline hazard", "Haz", "Baseline survival", "Su"), nomatch = 0)
  
  
  

  
  
  
  
  if (plot.type == 0) {
    stop("estimator must be 'Baseline hazard' (or 'Haz'), 'Baseline survival' (or 'Su')")
  }
  
  
  if(!(plot.type==3 || plot.type==4)){
  
    median=FALSE
  }
    
  

  
  if(!missing(conf.bands)){
    if (!is.logical(conf.bands) || length(conf.bands) != 1) {
      stop("'conf.bands' should be a logical value (TRUE or FALSE).")
    }
  }
  
  if((plot.type==3 || plot.type==4)){
  if(!missing(median)){
    if (!is.logical(median) || length(median) != 1) {
      stop("'median' should be a logical value (TRUE or FALSE).")
    }
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
  
  valid_transitions <- c("01","02","12")
  if(!missing(transition)){
    if (!(transition %in% valid_transitions)) {
      stop("'transition' should be one of: '01','02','12'.")
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
  
  
  
  
  ### HAZARD
  
  
  
  if(!missing(transition)){
    if(plot.type==1 || plot.type==2){
      if(conf.bands==TRUE){
        if(transition=="01"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 01"
          }
          no_conf <- NULL
          if(1 %in% partialH | 2 %in% partialH){
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 01 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
            matplot(x01, lam01[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(1 %in% partialH | 2 %in% partialH)) && dim(model$lam01)[2]==2 ){
            no_conf <- 1
            message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 01.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
            matplot(x01, lam01[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
          par(mfrow=c(1,1))
          lam01 <- hazard_and_confidence_bands_ID(x01,x$scale.weib[1],x$shape.weib[1],vcov[1:2, 1:2])
          matplot(x01, lam01[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        if(transition=="02"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 02"
          }
          no_conf <- NULL
          if(3 %in% partialH | 4 %in% partialH){
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 02 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
            matplot(x02, lam02[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(3 %in% partialH | 4 %in% partialH)) && dim(model$lam02)[2]==2 ){
            no_conf <- 1
            message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 02.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
            
            matplot(x02, lam02[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            par(mfrow=c(1,1))
            lam02 <- hazard_and_confidence_bands_ID(x02,x$scale.weib[2],x$shape.weib[2],vcov[3:4, 3:4])
            
            matplot(x02, lam02[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        
        if(transition=="12"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 12"
          }
          no_conf <- NULL
          if(5 %in% partialH | 6 %in% partialH){
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 12 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
            
            matplot(x12, lam12[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(5 %in% partialH | 6 %in% partialH)) && dim(model$lam12)[2]==2 ){
            no_conf <- 1
            lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
            message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 12.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            matplot(x12, lam12[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            par(mfrow=c(1,1))
            lam12 <- hazard_and_confidence_bands_ID(x12,x$scale.weib[3],x$shape.weib[3],vcov[5:6, 5:6])
            matplot(x12, lam12[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        
      }
      
      if(conf.bands==FALSE){
        if(transition=="01"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 01"
          }
          
          par(mfrow=c(1,1))
          lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
          matplot(x01, lam01[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = 1, cex = cex.legend)
        }
        
        
        if(transition=="02"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 02"
          }
          par(mfrow=c(1,1))
          lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
          matplot(x02, lam02[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = 1, cex = cex.legend)
        }
        
        
        if(transition=="12"){
          if(missing(Ylab)){
            Ylab <- "Baseline hazard function for transition 12"
          }
          par(mfrow=c(1,1))
          lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
          matplot(x12,lam12[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = 1, cex = cex.legend)
        }
      }
      
      
      
    }
  }
  
  
  
  
  
  if(missing(transition)){
    if(plot.type==1 || plot.type==2){
      if (conf.bands == TRUE) {
        par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
        if(missing(Ylab)){
          Ylab <- "Baseline hazard function"
        }
        if (1 %in% partialH | 2 %in% partialH) {
          lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
          message("Some parameters of the baseline hazard for transition 01 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x01, lam01[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(1 %in% partialH | 2 %in% partialH) && dim(model$lam01)[2] == 2) {
          lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
          message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 01.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x01, lam01[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
        } else {
          lam01 <- hazard_and_confidence_bands_ID(x01,x$scale.weib[1],x$shape.weib[1],vcov[1:2, 1:2])
          matplot(x01, lam01[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }
        
        if (3 %in% partialH | 4 %in% partialH) {
          lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
          message("Some parameters of the baseline hazard for transition 02 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x02, lam02[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(3 %in% partialH | 4 %in% partialH) && dim(model$lam02)[2] == 2) {
          lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
          message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 02.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x02, lam02[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else {
          lam02 <- hazard_and_confidence_bands_ID(x02,x$scale.weib[2],x$shape.weib[2],vcov[3:4, 3:4])
          matplot(x02, lam02[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }
        
        if (5 %in% partialH | 6 %in% partialH) {
          lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
          message("Some parameters of the baseline hazard for transition 12 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x12, lam12[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(5 %in% partialH | 6 %in% partialH) && dim(model$lam12)[2] == 2) {
          lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
          message("Numerical problem encountered in computing the confidence bands of the baseline hazard for transition 12.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x12, lam12[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else {
          lam12 <- hazard_and_confidence_bands_ID(x12,x$scale.weib[3],x$shape.weib[3],vcov[5:6, 5:6])
          matplot(x12, lam12[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline hazard", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }

      }
      
      if(conf.bands==FALSE){
        
        par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
        
        if(missing(Ylab)){
          Ylab <- "Baseline hazard function"
        }
        lam01 <- base_hazard_ID(x01,x$scale.weib[1],x$shape.weib[1])
        matplot(x01, lam01[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
        legend(pos.legend, legend = c("Baseline hazard"), 
               col = color, lty = 1, cex = cex.legend)
        
        
        lam02 <- base_hazard_ID(x02,x$scale.weib[2],x$shape.weib[2])
        
        matplot(x02, lam02[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
        legend(pos.legend, legend = c("Baseline hazard"), 
               col = color, lty = 1, cex = cex.legend)
        
        
        
        lam12 <- base_hazard_ID(x12,x$scale.weib[3],x$shape.weib[3])
        matplot(x12, lam12[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
        legend(pos.legend, legend = c("Baseline hazard"), 
               col = color, lty = 1, cex = cex.legend)
      }
      
      
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ### SURVIVAL
  
  if(!missing(transition)){
    if(plot.type==3 || plot.type==4){
      if(conf.bands==TRUE){
        if(transition=="01"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 01"
          }
          no_conf <- NULL
          if(1 %in% partialH | 2 %in% partialH){
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 01 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
            matplot(x01, surv01[,2], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(1 %in% partialH | 2 %in% partialH)) && dim(model$surv01)[2]==2 ){
            no_conf <- 1
            surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
            message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 01.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            matplot(x01, surv01[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            par(mfrow=c(1,1))
            surv01 <- survival_and_confidence_bands_ID(x01,x$scale.weib[1],x$shape.weib[1],vcov[1:2, 1:2])
            matplot(x01, surv01[,2:4], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        if(transition=="02"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 02"
          }
          no_conf <- NULL
          if(3 %in% partialH | 4 %in% partialH){
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 02 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
            matplot(x02, surv02[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(3 %in% partialH | 4 %in% partialH)) && dim(model$surv02)[2]==2 ){
            no_conf <- 1
            surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
            message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 02.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            matplot(x02, surv02[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            par(mfrow=c(1,1))
            surv02 <- survival_and_confidence_bands_ID(x02,x$scale.weib[2],x$shape.weib[2],vcov[3:4, 3:4])
            matplot(x02, surv02[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        
        if(transition=="12"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 12"
          }
          no_conf <- NULL
          if(5 %in% partialH | 6 %in% partialH){
            surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
            no_conf <- 1
            message("Some parameters of the baseline hazard for transition 12 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
            par(mfrow=c(1,1))
            matplot(x12, surv12[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          if((!(5 %in% partialH | 6 %in% partialH)) && dim(model$surv12)[2]==2 ){
            surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
            
            no_conf <- 1
            message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 12.")
            message("\nSuggestion: Try different initial values or increase the number of iterations.")
            matplot(x12, surv12[,2], col=color, type="l",  lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival"), 
                   col = color, lty = c(1), cex = cex.legend)
          }
          
          
          if(is.null(no_conf)){
            surv12 <- survival_and_confidence_bands_ID(x12,x$scale.weib[3],x$shape.weib[3],vcov[5:6, 5:6])
            
            par(mfrow=c(1,1))
            matplot(x12, surv12[,2:4], col=color, type="l", lty=c(1,2,2), xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
            legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                   col = color, lty = c(1, 2), cex = cex.legend)
          }
        }
        
        
      }
      
      if(conf.bands==FALSE){
        if(transition=="01"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 01"
          }
          
          par(mfrow=c(1,1))
          surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
          
          matplot(x01, surv01[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = 1, cex = cex.legend)
        }
        
        
        if(transition=="02"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 02"
          }
          par(mfrow=c(1,1))
          surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
          matplot(x02, surv02[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = 1, cex = cex.legend)
        }
        
        
        if(transition=="12"){
          if(missing(Ylab)){
            Ylab <- "Baseline survival function for transition 12"
          }
          par(mfrow=c(1,1))
          surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
          matplot(x12, surv12[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = 1, cex = cex.legend)
        }
      }
      
      
      
    }
  }
  
  
  
  
  
  if(missing(transition)){
    if(plot.type==3 || plot.type==4){
      if (conf.bands == TRUE) {
        par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
        if(missing(Ylab)){
          Ylab <- "Baseline survival function"
        }
        if (1 %in% partialH | 2 %in% partialH) {
          surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
          message("Some parameters of the baseline hazard for transition 01 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x01, surv01[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(1 %in% partialH | 2 %in% partialH) && dim(model$surv01)[2] == 2) {
          surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
          message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 01.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x01, surv01[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
        } else {
          surv01 <- survival_and_confidence_bands_ID(x01,x$scale.weib[1],x$shape.weib[1],vcov[1:2, 1:2])
          matplot(x01, surv01[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 01",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }
        
        if (3 %in% partialH | 4 %in% partialH) {
          surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
          message("Some parameters of the baseline hazard for transition 02 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x02, surv02[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(3 %in% partialH | 4 %in% partialH) && dim(model$surv02)[2] == 2) {
          surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
          message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 02.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x02, surv02[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
        } else {
          surv02 <- survival_and_confidence_bands_ID(x02,x$scale.weib[2],x$shape.weib[2],vcov[3:4, 3:4])
          matplot(x02, surv02[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 02",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }
        
        if (5 %in% partialH | 6 %in% partialH) {
          surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
          message("Some parameters of the baseline hazard for transition 12 were dropped from hessian in 'PartialH', confidence bands cannot be calculated in this case.")
          matplot(x12, surv12[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        } else if (!(5 %in% partialH | 6 %in% partialH) && dim(model$surv12)[2] == 2) {
          surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
          message("Numerical problem encountered in computing the confidence bands of the baseline survival for transition 12.")
          message("\nSuggestion: Try different initial values or increase the number of iterations.")
          matplot(x12, surv12[, 2], col=color, type="l", lty=1, xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
        } else {
          surv12 <- survival_and_confidence_bands_ID(x12,x$scale.weib[3],x$shape.weib[3],vcov[5:6, 5:6])
          matplot(x12, surv12[, 2:4], col=color, type="l", lty=c(1, 2, 2), xlab=Xlab, ylab=Ylab, main="Transition 12",lwd=lwd)
          legend(pos.legend, legend = c("Baseline survival", "Confidence bands"), 
                 col = color, lty = c(1, 2), cex = cex.legend)
        }

      }
      
      if(conf.bands==FALSE){
        
        par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
        
        if(missing(Ylab)){
          Ylab <- "Baseline survival function"
        }
        surv01 <- su_ID(x01,x$scale.weib[1],x$shape.weib[1])
        matplot(x01, surv01[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 01",lwd=lwd)
        legend(pos.legend, legend = c("Baseline survival"), 
               col = color, lty = 1, cex = cex.legend)
        
        
        
        surv02 <- su_ID(x02,x$scale.weib[2],x$shape.weib[2])
        matplot(x02, surv02[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 02",lwd=lwd)
        legend(pos.legend, legend = c("Baseline survival"), 
               col = color, lty = 1, cex = cex.legend)
        
        
        
        surv12 <- su_ID(x12,x$scale.weib[3],x$shape.weib[3])
        matplot(x12, surv12[,2], col=color, type="l", lty=1, xlab=Xlab,ylab=Ylab, main="Transition 12",lwd=lwd)
        legend(pos.legend, legend = c("Baseline survival"), 
               col = color, lty = 1, cex = cex.legend)

      }
      
      
      
    }
  }
  
  return(invisible())
}





