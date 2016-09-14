# Test statistic functions
Max_Quan_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, scale=FALSE) {
  # Computes maximum difference of quantiles 
  # Comments: OBC Using Linear Interpolation
  # of the ECDF Recall: ECDF range is [1/n,1]. This is probably not realistic 
  # Args: 
  # x: A vector of observations for a R.V. (must be numeric) 
  # y: Either (1) Another vectorof observations (two sample) 
  # (2) A quantile function such as qnorm -- must start with q (one sample)
  # interp: method of interpolation used. For more details, see ?quantile 
  # do.plot: Creates a plot illustrating the statistic 
  # Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
    }
    leny <- length(y)
    y <- sort(y)
    # x1 and y1 are the quantile values
    x1 <- seq(1/lenx, 1, 1/lenx)
    y1 <- seq(1/leny, 1, 1/leny)
    if (lenx == leny) {
      # In the case of equal sample sizes take max of abs. val
      z <- max(abs(y - x))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    } else if (lenx > leny) {
      # If there are unequal sample sizes, then interpolation is necessary
      q1 <- quantile(x, probs = x1, type = interp)
      q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
      z <- max(abs(q_inter(y1) - y))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    } else {
      # So length y>x
      q2 <- quantile(y, probs = y1, type = interp)
      q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
      z <- max(abs(q_inter(x1) - x))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    }
    return(z)
  }
  
  # One Sample
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  z <- max(abs(x - y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...)))  #Note: quantiles up for debate
  if (do.plot == TRUE) {
    plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
  }
  return(z)
}

Trap_Quan_Area_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, size = 0.25) {
  # Computes area based on trapezoid areas Args: x: A vector of observations for a
  # R.V. (must be numeric) y: Either (1) Another vector of observations (two sample)
  # (2) A quantile function such as qnorm (one sample) size: controls height of
  # trapezoids Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  
  # One Sample
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  # Quantiles
  x1 <- seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx)
  # Parameter which determines height of trapezoids
  delta <- size * (1/lenx)
  # Creating quantile function q1 <- quantile(x, probs = x1, type = interp) q_inter <-
  # approxfun(x1, q1, yleft = min(q1), yright = max(q1)) Length of horizontal segments
  x_p_d <- x + (y(x1 + delta, ...) - y(x1, ...))
  x_m_d <- x + (y(x1 - delta, ...) - y(x1, ...))
  y_p_d <- y(x1 + delta, ...)
  y_m_d <- y(x1 + delta, ...)
  # sum of trapezoid areas
  z <- sum((abs(x_p_d - y_p_d) + abs(x_m_d - y_m_d)) * delta)
  if (do.plot == TRUE) {
    plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
    
  }
  return(z)
}
Com_KS_Max_Quan_TS <- function(x, y, ...) {
  obc.stat <- Max_Quan_TS(x, y)
  ks.stat <- ks.res.simp(x, y)
  if (obc.stat > ks.stat) {
    z <- list(obc.stat, c("OBC"))
  } else if (ks.stat > obc.stat) {
    z <- list(ks.stat, c("KS"))
  } else z <- list(c(NULL), c("NEITHER"))
  
  return(z)
}
Com_KS_Max_Quan_TS_Simp <- function(x, y, ...) {
  obc.stat <- myts(x, y, ...)
  ks.stat <- ks.res.simp(x, y, ...)
  if (obc.stat > ks.stat) {
    z <- obc.stat
  } else if (ks.stat > obc.stat) {
    z <- ks.stat
  } else z <- 1000  #######FIX THIS LATER
  
  return(z)
}
Com_KS_Max_Quan_Range_TS <- function(x, y) {
  # Returns value of the statistic and whether or not it from OBC
  com <- c(x, y)
  range <- max(com) - min(com)
  obc.stat <- Max_Quan_TS(x, y)/range
  ks.stat <- ks.res(x, y)
  if (obc.stat > ks.stat) {
    z <- list(obc.stat, c("OBC"))
  } else if (ks.stat > obc.stat) {
    z <- list(ks.stat, c("KS"))
  } else z <- list(c(NULL), c("NEITHER"))
  
  return(z)
}
Com_KS_Max_Quan_Range_TS_Simp <- function(x, y) {
  # Returns value of the statistic and whether or not it from OBC
  com <- c(x, y)
  range <- max(com) - min(com)
  obc.stat <- Max_Quan_TS(x, y)/range
  ks.stat <- ks.res.simp(x, y)
  z <- max(obc.stat, ks.stat)
  return(z)
}

Density_Estimate <- function(x,interp=4){
  # Take in data (x)
  # Output: Function 
  # Outdated, see Kernel_CDF_Estimate
  lenx <- length(x)
  p_x <- seq(1/lenx, 1, 1/lenx)
  q1 <- quantile(x, probs = p_x, type = interp)
  dens_est <- approxfun(p_x, q1, yleft = min(q1), yright = max(q1))
}

Kernel_CDF_Estimate <- function(x, opt="cumsum"){
  if(opt=="integrate"){
    pdf <- density(x)
    f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
    cdf_fun <- function(lim) {
      cdf <- integrate(f, -Inf, lim, stop.on.error=FALSE)
    }
    return(cdf_fun)
  }
  if(opt=="cumsum"){
    pdf <- density(x)
    y <- cumsum(pdf$y)
    cdf <- pdf
    cdf$y <- y/max(y)
    cdf <- approxfun(cdf$x, cdf$y, yleft=0, yright=1)
    return(cdf)
  }
}

Kernel_CDF_Inverse <- function(x, opt="cumsum"){
  if(opt=="integrate"){
    pdf <- density(x)
    f <- approxfun(pdf$y, pdf$x, yleft=0, yright=0)
    inv_cdf_fun <- function(lim) {
      cdf <- integrate(f, -Inf, lim, stop.on.error=FALSE)
    }
    return(inv_cdf_fun)
  }
  if(opt=="cumsum"){
    pdf <- density(x)
    y <- cumsum(pdf$y)
    inv_cdf <- pdf
    inv_cdf$y <- y/max(y)
    inv_cdf <- approxfun(inv_cdf$y, inv_cdf$x, yleft=min(inv_cdf$x), yright=max(inv_cdf$x))
    return(inv_cdf)
  }
  
}

Kernel_Estimates <- function(x, opt="cumsum"){
  # Efficient combining of Kernel_CDF_Inverse and Kernel_CDF_Estimate
  if(opt=="cumsum"){
    pdf <- density(x)
    y <- cumsum(pdf$y)
    cdf <- pdf
    cdf$y <- y/max(y)
    inv_cdf <- approxfun(cdf$y, cdf$x, yleft=min(cdf$x), yright=max(cdf$x))
    cdf <- approxfun(cdf$x, cdf$y, yleft=0, yright=1)
    return(list("CDF" = cdf,"CDF Inverse"=inv_cdf))
  }
}

Delta_Calc <- function(x, y, lenx, leny){
  x_min <- min(abs(diff(x)))
  y_min <- min(abs(diff(y)))
  delta <- min(x_min, y_min)
}

Quad_OS_Area_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, opt=FALSE) {
  # Computes area based quadrilateral areas based on the order statistics
  # Calculates n boxes, where n is the sample size
  #
  # TODO: make work for unequal sample sizes
  #
  # Args: x: A vector of observations for a R.V. (must be numeric) 
  # y: Either (1) Another vector of observations (two sample)
  # Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  y <- sort(y)
  leny <- length(y)
  
  # Placeholder for Delta
  delta <- Delta_Calc(x,y)
  
  # Two Sample Test
  
  # Density estimate for x and y samples
  x_density <- Kernel_CDF_Estimate(x)
  y_density <- Kernel_CDF_Estimate(y)
  
  # Calculating Test Statistic
  x_minus_delta <- x - delta
  x_plus_delta <- x + delta
  
  #p_x_minus_delta <- x_density(x_minus_delta)$value
  p_x_minus_delta <- sapply(x_minus_delta, function(x) x_density(x))
  
  #p_x_plus_delta <- x_density(x_plus_delta)$value
  p_x_plus_delta <- sapply(x_plus_delta, function(x) x_density(x))
  
  y_minus_delta <- y - delta
  y_plus_delta <- y + delta
  
  #p_y_minus_delta <- y_density(y_minus_delta)$value
  p_y_minus_delta <- sapply(y_minus_delta, function(y) y_density(y))
  
  #p_y_plus_delta <- y_density(y_plus_delta)$values
  p_y_plus_delta <- sapply(y_plus_delta, function(y) y_density(y))
  
  x1 <- as.list(as.data.frame(rbind(x_minus_delta, p_x_minus_delta)))
  x2 <- as.list(as.data.frame(rbind(x_plus_delta, p_x_plus_delta)))
  y1 <- as.list(as.data.frame(rbind(y_minus_delta, p_y_minus_delta)))
  y2 <- as.list(as.data.frame(rbind(y_plus_delta, p_y_plus_delta)))
  
  coords <- rbind(x1,x2,y2,y1)
  coords <- as.list(as.data.frame(coords))

  # Draw Plot Option
  if(do.plot==TRUE){
    Plot_Quad_Areas(x=x, y=y, lenx=lenx, coords=coords, xdens=x_density, ydens=y_density)
  }
  # Diagnostic Option. Can be deleted later.
  if(opt=="coords"){
    return(coords)
  }
  # Calculating the value of the statistic
  areas <- lapply(coords, function(x) quad.area(x[[1]],x[[2]],x[[3]],x[[4]]))
  total_area <- sum(unlist(areas))
}

Quad_Quan_Area_TS <- function(x, y, ..., interp = 4, 
                              do.plot = FALSE, 
                              opt=FALSE,
                              maxval=FALSE, 
                              div=4) {
  # Computes area based quadrilateral areas based on the order statistics
  # Calculates n boxes, where n is the sample size
  #
  # TODO: make work for unequal sample sizes
  #
  # Args: x: A vector of observations for a R.V. (must be numeric) 
  # y: Either (1) Another vector of observations (two sample)
  # Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  y <- sort(y)
  leny <- length(y)
  
  # Placeholder for Delta
  #delta <- Delta_Calc(x,y)
  
  # Two Sample Test
  
  # Quantile function estimate for x and y samples
  x_functions <- Kernel_Estimates(x)
  y_functions <- Kernel_Estimates(y)

  x_quantile_fun <- x_functions[["CDF Inverse"]]
  y_quantile_fun <- y_functions[["CDF Inverse"]]
  
  x_density <- x_functions[["CDF"]]
  y_density <- y_functions[["CDF"]]
  
  
  # Finding values of the quantile function desired
  # Should be leq to n
  n <- round(min(lenx, leny) / div)
  p <- seq(1/(n+1), n/(n+1), length.out=n) 
    
  # Start from the y axis (0 < p < 1)
  x_quantiles <- x_quantile_fun(p)
  y_quantiles <- y_quantile_fun(p)

  # New Delta Calc
  delta <- min(diff(x_quantiles)/2, diff(y_quantiles)/2)
  
  # Now that we are on data axis, add plus or minus delta
  x_minus_delta <- x_quantiles - delta
  x_plus_delta <- x_quantiles + delta
  y_minus_delta <- y_quantiles - delta
  y_plus_delta <- y_quantiles + delta
  
  # Reapply density estimates
  p_x_minus_delta <- sapply(x_minus_delta, function(x) x_density(x))
  p_x_plus_delta <- sapply(x_plus_delta, function(x) x_density(x))
  p_y_minus_delta <- sapply(y_minus_delta, function(x) y_density(x))
  p_y_plus_delta <- sapply(y_plus_delta, function(x) y_density(x))
  
  # getting everything as coordinates
  # (x,y)
  x1 <- as.list(as.data.frame(rbind(x_minus_delta, p_x_minus_delta)))
  x2 <- as.list(as.data.frame(rbind(x_plus_delta, p_x_plus_delta)))
  y1 <- as.list(as.data.frame(rbind(y_minus_delta, p_y_minus_delta)))
  y2 <- as.list(as.data.frame(rbind(y_plus_delta, p_y_plus_delta)))
  
  coords <- rbind(x1,x2,y2,y1)
  # (x - delta, f(x - delta))
  # (x + delta, f(x + delta))
  # (y - delta, f(y - delta))
  # (y + delta, f(y + delta)) 
  coords <- as.list(as.data.frame(coords))

  # Diagnostic Option. Can be deleted later.
  if(opt=="coords"){
    return(coords)
  }
  # Calculating the area values
  areas <- lapply(coords, function(x) quad.area(x[[1]],x[[2]],x[[3]],x[[4]]))
  
  # Draw Plot Option
  if(do.plot==TRUE & maxval==FALSE){
    Plot_Quad_Areas(x=x,y=y,lenx=lenx,coords=coords, xdens=x_density, ydens=y_density)
  }
  if(do.plot==TRUE & maxval==TRUE){
    Plot_Quad_Areas(x=x,y=y,lenx=lenx,coords=coords, xdens=x_density, ydens=y_density, areas=areas)
  }
  if(maxval==TRUE){
    return(max(unlist(areas)))
  }  
  # Calculating total area
  total_area <- sum(unlist(areas))
}

