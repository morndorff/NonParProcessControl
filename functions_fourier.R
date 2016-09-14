# Developing One Sample Test based on LaRiccia

fourier_test <- function(x, y, ..., doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  
  # Is equivalent to CvM
  
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=lenx)
  for(j in 1:lenx){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  
  C_n <- a_jn^2 * (1 / ((1:lenx) * pi)^2 )
  if(doplot){
    plot(C_n, type="l", xlab="J Index", ylab="Component Value")
    plot(cumsum(C_n), type="l", xlab="J Index", ylab="Cum Sum")
  }
  C_n <- lenx * sum(C_n)
  C_n
}


CVM_Fourier <- function(x, y, ..., doplot=FALSE, type="durbin", num_coef=500){
  # Args: x - vector
  # y : pnorm or similar
  
  # Is equivalent to CvM
  
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)

  # OPTION DURBIN IS CORRECT
  if(type=="durbin"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    for(j in 1:num_coef){
      zn[j] <- (sqrt(2/lenx)) * sum(cos(j*pi*x_unit))
      Wn[j] <- zn[j]^2 / ((j^2)*pi^2)
    }
    STAT <- sum(Wn)
    return(STAT)
  }
  
  if(type=="fanortho"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    #print(x_unit)
    for(j in 1:num_coef){
      theta1 <- mean(cos(2*j*pi*x_unit))
      theta2 <- mean(sin(2*j*pi*x_unit)) # (sqrt(2/lenx)) * sum(sin(2*j*pi*x_unit))
      zn[j] <- theta1^2 + theta2^2
      #print(zn)
      Wn[j] <-  zn[j] / (j^2)
      
    }
    #print(Wn)
    STAT <- sum(Wn) * lenx / (2*pi^2)
    return(STAT)
  }
  # Fan Attempt2
  if(type=="fanortho2"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    #print(x_unit)
    for(j in 1:num_coef){
      theta1 <- 1/(sqrt(2)*lenx) * sum(cos(2*j*pi*x_unit))
      theta2 <- 1/(sqrt(2)*lenx) * sum(sin(2*j*pi*x_unit))
      zn[j] <- theta1^2 + theta2^2
      Wn[j] <- zn[j] / (j^2)
    }
    #print(Wn)
    STAT <- sum(Wn) * lenx / (2*pi^2)
    return(STAT)
  }
  # How I would do this
  if(type=="mark"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    #print(x_unit)
    for(j in 1:num_coef){
      theta1 <- (1 / (j * pi * sqrt(lenx))) * sum(cos(2*j*pi*x_unit))
      theta2 <- (1 / (j * pi * sqrt(lenx))) * sum(sin(2*j*pi*x_unit))
      zn[j] <- theta1^2 + theta2^2
      Wn[j] <- zn[j] 
    }
    #print(Wn)
    STAT <- sum(Wn) + mean(x_unit)#* lenx / (2*pi^2)
    return(STAT)
  }
  if(type=="fft"){
    fx <- fft(x_unit)
    
    cos_p <- Re(fx)
    sin_p <- Im(fx)
    cos_p <- cos_p[-1]
    print(sin_p)
    num <- vector(length=length(cos_p))
    for(j in seq_along(cos_p)){
      num[j] <- (sin_p[j]^2 + cos_p[j]^2)/j^2
    }
    STAT <- sum(num) 
    return(STAT)
  }
  
  if(type=="chicken"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    
    #x_unit = F0(x)
    F_x <- ecdf(x_unit)
    y <- F_x(x_unit)
    newx <- y - x_unit
    print(newx)
    for(j in 1:num_coef){
      theta1 <- mean(cos(2*j*pi*newx))
      theta2 <- mean(sin(2*j*pi*newx)) # (sqrt(2/lenx)) * sum(sin(2*j*pi*x_unit))
      zn[j] <- theta1^2 + theta2^2
      #print(zn)
      Wn[j] <-  zn[j] / (j^2)
    }
    #print(Wn)
    STAT <- sum(Wn) * lenx / (2*pi^2)
    return(STAT)
  }
  
  if(type=="fan"){
    zn <- vector(length=num_coef)
    Wn <- vector(length=num_coef)
    #print(x_unit)
    for(j in 1:num_coef){
      theta1 <-  (sqrt(2/lenx)) * sum(cos(2*j*pi*x_unit))
      theta2 <- (sqrt(2/lenx)) * sum(sin(2*j*pi*x_unit))
      zn[j] <- theta1^2 + theta2^2
      #print(zn)
      Wn[j] <-  zn[j] / (j^2)
      
    }
    #print(Wn)
    STAT <- sum(Wn) * lenx / (2*pi^2)
    return(STAT)
  }
  
  
  if(type=="chicken2"){
  F_x <- ecdf(x_unit)
  y <- F_x(x_unit)
  newx <- y - x_unit
  #theta <- theta[-1]
  #jcos <- Re(theta) / lenx
  #jsin <- Im(theta) / lenx
  track <- vector(mode="numeric", length=length(jcos))
  for(i in 1:length(jcos)){
    track[i] <- jcos[i]^2 + jsin[i]^2
    track[i] <- track[i] / i^2
  }
  STAT <- (lenx / 2*pi^2) * sum(track)
  return(STAT)
  }
}






smooth_test <- function(x, y, ..., M=5, doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=M)
  for(j in 1:M){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  T_nm <- lenx * sum(a_jn^2)
  T_nm
}

smooth_test2 <- function(x, y, ..., lambda=5, doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  # lambda: smoothing parameter > 0
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=lenx)
  for(j in 1:lenx){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  seq_n <- 1:lenx
  S_lambda <- lenx * sum(a_jn^2 / (1 + (lambda * seq_n^2))^2)
  S_lambda
}

fft_test <- function(x,y, ...){
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  fft(x_unit)
}

Fan_Fourier_Test <- function(x, y, ..., m=2){
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  
  # Creating Sub Intervals
  nbin <- 2^(floor(log2(lenx - 1))) # Recommended by Fan
  # print(nbin)
  # bp <- seq(0, 1, length.out=nbin) # Bins are equally spaced
  #print(x_unit)
  # int <- findInterval(x_unit, bp) # Binning
  # counts <- as.numeric(table(int))   # Gives counts in each interval
  # print(counts)
  
  br <- seq(0, 1, length.out=nbin + 1)
  counts <- hist(x_unit, br, plot=F)$counts
  print(counts)
  
  Y_j <- 2 * (sqrt(counts)- sqrt(lenx/nbin)) # Square root transform
  # print(" Bin Transform")
  #print(Y_j)
  # Now is Y ~ N (theta, Inbin)
  # sum(Y_j)
  #print(Y_j)
  # Now we take orthogonal transform
  X_lambda <- fft(Y_j)
  #print("fft")
  # print(X_lambda)
  a <- Re(X_lambda[2:length(X_lambda)]) / lenx # Discard scaling coefficnet
  b <- Im(X_lambda[2:length(X_lambda)]) / lenx
  
  #X_lambda <- (abs(X_lambda)/ lenx)^2
  #print("squared")
  #print(X_lambda)
  #print(X_lambda)
  STAT <- sum(a[1:m]^2, b[1:m]^2)
  STAT
}
