wave.den <- function (x, y, ..., doplot=F, wf="haar", scale=FALSE)  #n=2^5
{
  # Inputs:
  # Y must either be numeric or "p"dist
  F.x <- ecdf(x)
  
  # Two Sample 
  if(is.numeric(y)){
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      F.x <- ecdf(x)
      # return(newdata)
    }
    F.y <- ecdf(y)
    ml <- min(length(x),length(y)) 
    n <- 2^floor(log2(ml))
    z <- seq(range(x, y)[1], range(x, y)[2], length.out=n)
    F.dwt <- dwt(F.x(z) - F.y(z), wf=wf, n.levels=log(n, 2))
  } else{
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
    
    n <- 2^floor(log2(length(x)))
    z <- seq(min(x),max(x),length=n)
    F.dwt <- dwt(F.x(z) - y(z,...), wf=wf, n.levels=log(n,2))
  }
  oc <- unlist(F.dwt)
  oc <- max(abs(oc))
  #test_ks <- max(abs(F.x(z)-F.y(z)))
  
  if(doplot)
  {
    plot(z, F.x(z), type="l", 
         ylim=range(1.1, F.x(z), F.y(z), F.x(z) - F.y(z)))
    lines(z, F.y(z), col=3)
    lines(z, F.x(z) - F.y(z), col=4)
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    #segments(z[3], 0, z[3], ks, lwd=2)
    title(paste("oc = ", round(oc, 2),sep=""))
  }
  oc
}


wave.energy <- function (x, y, ...,
                         #n=2^5, 
                         doplot=F, 
                         opt="sum", 
                         wf="haar",
                         square=FALSE,
                         norm=TRUE,
                         scale=FALSE)
{
  # Args:
  # y must be numeric of a "p" distribution function
  if(opt !="max" & opt != "sum") stop ("invalid option value")
  # Get cdfs:
  F.x <- ecdf(x)
  # Two Sample
  if(is.numeric(y)){
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      F.x <- ecdf(x)
    }
    F.y <- ecdf(y)
    
    ml <- min(length(x),length(y))
    n <- 2^floor(log2(ml))
    z <- seq(range(x, y)[1], range(x, y)[2], length=n)
    
    F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
    F.y.dwt <- dwt(F.y(z), wf=wf, n.levels=log(n, 2))
    x.dwt <- unlist(F.x.dwt)
    y.dwt <- unlist(F.y.dwt)
  }else{
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
    
    n <- 2^floor(log2(length(x)))
    z <- seq(min(x), max(x), length=n)
    F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
    F.y.dwt <- dwt(y(z,...), wf=wf, n.levels=log(n, 2))
    x.dwt <- unlist(F.x.dwt)
    y.dwt <- unlist(F.y.dwt)
  }
  
  if(square==FALSE){
    x.dwt <- sort(abs(x.dwt), decreasing = TRUE)
    y.dwt <- sort(abs(y.dwt), decreasing = TRUE)
  }
  # What if we square the coefficients?
  if(square){
    x.dwt <- sort(x.dwt^2, decreasing=TRUE)
    y.dwt <- sort(y.dwt^2, decreasing=TRUE)
  }
  
  xx <- cumsum(x.dwt)
  yy <- cumsum(y.dwt)
  max_x <- tail(xx, 1)
  max_y <- tail(yy, 1)
  # return(max(yy))
  if(norm==FALSE){
    if(max_x > max_y){
      xx <- xx / max_x
      yy <- yy / max_x
    } else{
      xx <- xx / max_y
      yy <- yy / max_y
    }
  }
  if(norm==TRUE){  
    xx <- xx / max_x
    yy <- yy / max_y
  }
  if(doplot){
    plot(xx, type="l", ylim=range(xx, yy, xx - yy))
    lines(yy, col=2)
    lines(xx - yy, col=3)
    abline(h=0, lty=3)
  }
  
  # Returning test satistic values
  if(opt=="max"){
    ts <- max(abs(xx-yy))
    names(ts) <- "max"
    return(ts)
  }
  if(opt=="sum"){
    ts <- sum(abs(xx-yy))
    names(ts) <- "sum"
    return(ts)
  }
}


wave.bec <- function(x,y, ..., interp = 4, doplot=F, wf="haar", reduce=2, scale=FALSE)
{
  library(waveslim)
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
    }
    x <- sort(x)
    leny <- length(y)
    y <- sort(y)

    num_quan <- min(2^floor(log(lenx/reduce,2)), 2^floor(log(leny/reduce,2)))
    prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)

    # Because # of quantiles is < data points, need to interpolate
    qx <- quantile(x, probs = prob, type = interp)
    qy <- quantile(y, probs = prob, type = interp)
    q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
    q_inter_y <- approxfun(prob, qy, yleft = min(qy), yright = max(qy))
    quan_dif <- q_inter_y(prob) - q_inter_x(prob)
    
    test_wave <- dwt(quan_dif,wf="haar")
    w_coef <- unlist(test_wave)
    # Take abs. value of wavelet coefficients, sort in decreasing order
    w_coef <- sort(abs(w_coef), decreasing=TRUE) 
    # Sum of them
    sum_coef <- sum(w_coef)
    # Get number of coefficients to keep, based on 90% thresholding
    num_coef <- which(cumsum(abs(w_coef))/sum_coef<=.9)
    th_w_coef<- w_coef[num_coef]
    STAT <- sum(th_w_coef^2)
    return(STAT)
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
  # Assuring diadic lengths
  #nx_p2 <- 2^floor(log(lenx,2))
  
  num_quan <- 2^floor(log(lenx/reduce,2))
  prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
    
  # Estimated Quantiles:
  qx <- quantile(x, probs = prob, type = interp)
  q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
  x_quans <- q_inter_x(prob)
  # True Quantiles
  true_quan <- y(prob, ...)
  quan_dif <- x_quans - true_quan
  
  test_wave <- dwt(quan_dif,wf="haar", n.levels=log(num_quan,2))
  w_coef <- unlist(test_wave)
  # Take abs. value of wavelet coefficients, sort in decreasing order
  w_coef <- sort(abs(w_coef), decreasing=TRUE) 
  # Sum of them
  sum_coef <- sum(w_coef)
  # Get number of coefficients to keep, based on 90% thresholding
  num_coef <- min(1, which(cumsum(abs(w_coef))/sum_coef<=.9))
  th_w_coef<- w_coef[num_coef]
  STAT <- sum(th_w_coef^2)
  return(STAT)
}

wave.energy2 <- function (x, y, ...,
                         n=2^4, 
                         doplot=F, 
                         opt="sum", 
                         wf="haar",
                         square=FALSE,
                         norm=TRUE,
                         scale=FALSE)
{
  # Args:
  # y must be numeric of a "p" distribution function
  if(opt !="max" & opt != "sum") stop ("invalid option value")
  # Get cdfs:
  F.x <- ecdf(x)
  # Two Sample
  if(is.numeric(y)){
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      F.x <- ecdf(x)
    }
    F.y <- ecdf(y)
    # if(n > min(length(x),length(y))) n <- 2^floor(log(min(length(x),length(y)),2))
    
    # ml <- min(length(x),length(y))
    # n <- 2^floor(log2(ml))
    z <- seq(range(x, y)[1], range(x, y)[2], length=n)
    
    F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
    F.y.dwt <- dwt(F.y(z), wf=wf, n.levels=log(n, 2))
    x.dwt <- unlist(F.x.dwt)
    y.dwt <- unlist(F.y.dwt)
  }else{
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
    # if(n > length(x)) n <- 2^floor(log(length(x),2))
    z <- seq(min(x), max(x), length=n)
    F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
    F.y.dwt <- dwt(y(z,...), wf=wf, n.levels=log(n, 2))
    x.dwt <- unlist(F.x.dwt)
    y.dwt <- unlist(F.y.dwt)
  }
  
  if(square==FALSE){
    x.dwt <- sort(abs(x.dwt), decreasing = TRUE)
    y.dwt <- sort(abs(y.dwt), decreasing = TRUE)
  }
  # What if we square the coefficients?
  if(square){
    x.dwt <- sort(x.dwt^2, decreasing=TRUE)
    y.dwt <- sort(y.dwt^2, decreasing=TRUE)
  }
  
  xx <- cumsum(x.dwt)
  yy <- cumsum(y.dwt)
  max_x <- tail(xx, 1)
  max_y <- tail(yy, 1)
  # return(max(yy))
  if(norm==FALSE){
    if(max_x > max_y){
      xx <- xx / max_x
      yy <- yy / max_x
    } else{
      xx <- xx / max_y
      yy <- yy / max_y
    }
  }
  if(norm==TRUE){  
    xx <- xx / max_x
    yy <- yy / max_y
  }
  if(doplot){
    plot(xx, type="l", ylim=range(xx, yy, xx - yy))
    lines(yy, col=2)
    lines(xx - yy, col=3)
    abline(h=0, lty=3)
  }
  
  # Returning test satistic values
  if(opt=="max"){
    ts <- max(abs(xx-yy))
    names(ts) <- "max"
    return(ts)
  }
  if(opt=="sum"){
    ts <- sum(abs(xx-yy))
    names(ts) <- "sum"
    return(ts)
  }
}

wave.den2 <- function (x, y, ..., doplot=F, wf="haar", n=2^5, scale=FALSE)
{
  # Inputs:
  # Y must either be numeric or "p"dist
  one <- FALSE; two <- FALSE
  F.x <- ecdf(x)
  # Two Sample 
  if(is.numeric(y)){
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      F.x <- ecdf(x)
    }
    if(n > min(length(x),length(y))) n <- 2^floor(log(min(length(x),length(y)),2))
    F.y <- ecdf(y)
    ml <- min(length(x),length(y)) 
    # n <- 2^floor(log2(ml))
    z <- seq(range(x, y)[1], range(x, y)[2], length.out=n)
    F.dwt <- dwt(F.x(z) - F.y(z), wf=wf, n.levels=log(n, 2))
    two <- TRUE
  } else{
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
    if(n > length(x)) n <- 2^floor(log(length(x),2))
    z <- seq(min(x),max(x),length=n)
    F.dwt <- dwt(F.x(z) - y(z,...), wf=wf, n.levels=log(n,2))
    one <- TRUE
  }
  oc <- unlist(F.dwt)
  oc <- max(abs(oc))
  
  if(doplot & two){
    plot(z, F.x(z), type="l", 
         ylim=range(1.1, F.x(z), F.y(z), F.x(z) - F.y(z)))
    lines(z, F.y(z), col=3)
    lines(z, F.x(z) - F.y(z), col=4)
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    #segments(z[3], 0, z[3], ks, lwd=2)
    title(paste("oc = ", round(oc, 2),sep=""))
  }
  if(doplot & one){
    plot(z,F.x(z),type="l",
         ylim=range(1.1, F.x(z), y(z), F.x(z) - y(z)))
    lines(z, y(z,...), col=3)
    lines(z, F.x(z) - y (z), col=4)
    abline(h=0, lty=3)
    abline(h=0, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    title(paste("oc = ", round(oc, 2),sep=""))
  }

  #test_ks <- max(abs(F.x(z)-F.y(z)))
  oc

}

wave.bec2 <- function(x,y, ..., interp = 4, doplot=F, wf="haar", n=2^5, scale=FALSE)
{
  library(waveslim)
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      x <- sort(x)
    }
    leny <- length(y)
    y <- sort(y)
    if(n > min(length(x),length(y))) n <- 2^floor(log(min(length(x),length(y)),2))
    num_quan <- n
    # num_quan <- min(2^floor(log(lenx/reduce,2)), 2^floor(log(leny/reduce,2)))
    prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
    
    # Because # of quantiles is < data points, need to interpolate
    qx <- quantile(x, probs = prob, type = interp)
    qy <- quantile(y, probs = prob, type = interp)
    q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
    q_inter_y <- approxfun(prob, qy, yleft = min(qy), yright = max(qy))
    quan_dif <- q_inter_y(prob) - q_inter_x(prob)
    
    test_wave <- dwt(quan_dif,wf="haar")
    w_coef <- unlist(test_wave)
    # Take abs. value of wavelet coefficients, sort in decreasing order
    w_coef <- sort(abs(w_coef), decreasing=TRUE) 
    # Sum of them
    sum_coef <- sum(w_coef)
    # Get number of coefficients to keep, based on 90% thresholding
    num_coef <- which(cumsum(abs(w_coef))/sum_coef<=.9)
    th_w_coef<- w_coef[num_coef]
    STAT <- sum(th_w_coef^2)
    return(STAT)
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
  
  # Assuring diadic lengths
  #nx_p2 <- 2^floor(log(lenx,2))
  if(n > length(x)) n <- 2^floor(log(length(x),2))
  #num_quan <- 2^floor(log(lenx/reduce,2))
  
  num_quan <- n
  prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
  
  # Estimated Quantiles:
  qx <- quantile(x, probs = prob, type = interp)
  q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
  x_quans <- q_inter_x(prob)
  # True Quantiles
  true_quan <- y(prob, ...)
  quan_dif <- x_quans - true_quan
  
  test_wave <- dwt(quan_dif,wf="haar", n.levels=log(num_quan,2))
  w_coef <- unlist(test_wave)
  # Take abs. value of wavelet coefficients, sort in decreasing order
  w_coef <- sort(abs(w_coef), decreasing=TRUE) 
  # Sum of them
  sum_coef <- sum(w_coef)
  # Get number of coefficients to keep, based on 90% thresholding
  num_coef <- min(1, which(cumsum(abs(w_coef))/sum_coef<=.9))
  
  th_w_coef<- w_coef[num_coef]
  STAT <- sum(th_w_coef^2)
  return(STAT)
}

wave.mark <- function(x, y, ..., doplot=F, wf="haar", scale=FALSE)  #n=2^5
{
  # Inputs:
  # Y must either be numeric or "p"dist
  F.x <- ecdf(x)
  
  # Two Sample 
  if(is.numeric(y)){
    if(scale){
      newdata <- scale_two_sample(x,y)
      x <- newdata[[1]]
      y <- newdata[[2]]
      F.x <- ecdf(x)
      # return(newdata)
    }
    F.y <- ecdf(y)
    ml <- min(length(x),length(y)) 
    n <- 2^floor(log2(ml))
    z <- seq(range(x, y)[1], range(x, y)[2], length.out=n)
    F.dwt <- dwt(F.x(z) - F.y(z), wf=wf, n.levels=log(n, 2))
  } else{
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
    
    n <- 2^floor(log2(length(x)))
    z <- seq(min(x),max(x),length=n)
    F.dwt <- dwt(F.x(z) - y(z,...), wf=wf, n.levels=log(n,2))
    # Begin thresholding
    
  }
  oc <- unlist(F.dwt)
  oc <- max(abs(oc))
  #test_ks <- max(abs(F.x(z)-F.y(z)))
  
  if(doplot)
  {
    plot(z, F.x(z), type="l", 
         ylim=range(1.1, F.x(z), F.y(z), F.x(z) - F.y(z)))
    lines(z, F.y(z), col=3)
    lines(z, F.x(z) - F.y(z), col=4)
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    #segments(z[3], 0, z[3], ks, lwd=2)
    title(paste("oc = ", round(oc, 2),sep=""))
  }
  oc
}