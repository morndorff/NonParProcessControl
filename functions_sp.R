# SP Plot Functions
SP_Threshold_Selection <- function(x, y, ...,
                                grid_fine=15, dyadic_after=TRUE, 
                                peek_scaling=TRUE, plotEstFun=F, doplot=F,
                                debug=F, diag=F, plotThresh=F, includescale=F,
                                minthreshlevel=2, wf="haar"){
  # includescale: threshold lowest level scaling coefficent?
  # threshscale: Do we want to threshold the detail functions?
  # grid_fine: How fine the cross validation grid is
  # dyadic_after: Making the lambda value smaller or larger 
  # peek_scaling: Effects our scaling. Do we look at all the data when we scale,
  # or do we just look at the in sample data
  # If true, we do, if false we don't

  # Only works with dyadic sample size  
  require(wavethresh)
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  
  # One Sample -----
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  
  # Calculating Plot Differences -----
  
  x <- sort(x)
  lenx <- length(x)
  i <- 1:lenx
  sp <- (2/pi) * asin(sqrt((i-.5)/lenx))
  si <- (2/pi) * asin(sqrt(y(x, ...))) # Important part
  
  dif <- si-sp
  
  # Split the Samples -------
  if(lenx %% 2==0){
    #print("its even!")
    
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    si_even <- dif[ind_even]
    len_even <- lenx / 2
    si_odd <- dif[ind_odd]
    len_odd <- len_even
  } else{
    #print("its odd!")
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    
    si_even <- dif[ind_even]
    len_even <- (lenx - 1) / 2
    
    si_odd <- dif[ind_odd]
    len_odd <- len_even + 1
  } 
  # Getting Odd and Even Lambda -----
  Odd_OOS <- Find_Lambda_SP(in_sample_dif = si_even, out_of_sample_dif = si_odd,
                         grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot,
                         insamplefirst=F, plotThresh=plotThresh, 
                         includescale=includescale, minthreshlevel=minthreshlevel,
                         wf=wf)
  if(plotThresh) return(Odd_OOS)
  Odd_MSE <- Odd_OOS[["MSE Vector"]]
  Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]

  Even_OOS <- Find_Lambda_SP(in_sample_dif = si_odd, out_of_sample_dif =si_even,
                          grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot,
                          insamplefirst=T, plotThresh=plotThresh, 
                          includescale=includescale, minthreshlevel=minthreshlevel,
                          wf=wf)
  Even_MSE <- Even_OOS[["MSE Vector"]]
  Even_Lambda <- Even_OOS[["Optimal Lambda"]]

  # Choosing Lambda----------
  
  Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
  # This Chosen Lambda value needs to be adjusted to 
  # account for the fact that the even and odd samples
  # are of sample size n/2
  Average_Relative_Position <- mean(c(Even_OOS[["Position"]], Odd_OOS[["Position"]]))

  if(dyadic_after){
    # If the sample is going to be adjusted down or up to a
    # dyadic power of 2 afterwards, we should account for that
    Adjusted_Length <- 2^floor(log2(lenx))
    Adjusted_Lambda <- (1/ sqrt((1 - (log(2)/log(Adjusted_Length))))) * Chosen_Lambda
  } else{
    Adjusted_Lambda <- (1/ sqrt((1 - (log(2)/log(lenx))))) * Chosen_Lambda
  }
  if(diag){
    return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
                "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
                "Chosen Threshold (Lambda)"= Chosen_Lambda,
                "Sample Size Adjusted Threshold" = Adjusted_Lambda
                ))
  }
  return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
              "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
              "Chosen Threshold (Lambda)"= Chosen_Lambda,
              "Sample Size Adjusted Threshold"=Adjusted_Lambda,
              "Relative Position Chosen"=Average_Relative_Position))
}

Find_Lambda_SP <- function(in_sample_dif, out_of_sample_dif, 
                        grid_fine, plotEstFun, doplot, insamplefirst, plotThresh, includescale,
                        minthreshlevel=2, wf="haar"){
  # internal function for cross validation. Not to be used 
  
  # Input: A 'in sample' x and y and an out of sample x and y
  # Output: A list containing the MSE's used as well as the optimal lambda value found

  # minthreshlevel - At what level > 0  do we want to threshold the detail coefficients
  # includescale - Include the scaling coefficient in the thresholding? Only an option if minthreshlevel=0 

  if(includescale && minthreshlevel!=0) stop("Don't want to threshold scale function and only top level details")
  # Decompose input function
  if(wf=="la8") model_wd <- wd(in_sample_dif, filter.number=8, family="DaubLeAsymm") #la8
  if(wf=="haar") model_wd <- wd(in_sample_dif, filter.number=1, family="DaubExPhase") #haar
  
  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=minthreshlevel:(Largest_Level-1)))
  if(includescale){
    scaling <-  accessC(model_wd, level=0)
  }else{
    scaling <- NULL
  }
  Largest_Coefficient <- max(abs(c(details, scaling))) 
  Smallest_Coefficient <- min(abs(c(details, scaling)))
  
  lambda_grid <- seq(Smallest_Coefficient, 
                     Largest_Coefficient, length.out=grid_fine)
  
  
  # Making sure last threshold thresholds everything
  # lambda_grid <- append(lambda_grid, Largest_Coefficient + .1*sd(lambda_grid))
  
  # For values in the lambda grid,
  # Perform thresholding, reconstruction, and MSE comparison
  MSE <- vector(length=grid_fine)
  #Percent_Thresholded <- vector(length=grid_fine)
  Thresh_list <- vector(mode="list", length=grid_fine)
  for(i in seq_along(lambda_grid)){
    # Thresholding
    thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                             value=lambda_grid[i], 
                             levels=minthreshlevel:(nlevelsWT(model_wd)-1) )#, verbose=TRUE)
    # Percentage of Eligible Detail Coefficients Discarded 
    # details <- unlist(testfun(thresh_test, level=minthreshlevel:(Largest_Level-1)))
    # Percent_Thresholded[i] <- sum(details == 0) /length(details)
    
    if(includescale){
      Lowest_C <- abs(accessC(thresh_test, level=0))
      if(Lowest_C <= lambda_grid[i]){
        thresh_test <- putC(thresh_test, level=0, v=0)
      }
    }
  
    # Reconstruction
    reconstructed_fun <- wr(thresh_test)
    # Here, we estimate the odd from the even, and vice versa
    
    # First, we need to do some endpoint correction
    # Then, we take our estimate as the midpoint between entries
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    if(insamplefirst){
      Dif_First <- reconstructed_fun[1] / 2 # Take midpoint  
      estimated_fun <- ma(reconstructed_fun)
      estimated_fun <- c(Dif_First, head(estimated_fun, -1))
    }else{
      #Interpolate maximum
      Dif_End <- tail(reconstructed_fun , 1) / 2
      estimated_fun <- ma(reconstructed_fun)
      estimated_fun <- c(head(estimated_fun, -1), Dif_End)
    }
    
    # calculating MSE
    errors <- estimated_fun - out_of_sample_dif
    MSE[i] <- sum(errors^2)
    if(plotThresh) Thresh_list[[i]] <- (list("Lambda Value"=paste(i, "of", length(lambda_grid)), 
                               "Thresholded Coefficients"=thresh_test, 
                               "Reconstructed Function"=estimated_fun, 
                               "Unthresholded Wavelet"=model_wd))
    
    if(plotEstFun){
      #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
      par(mfrow=c(2,1))
      plot(in_sample_dif, type="l",
           main=paste("Lambda Position=", round(i/length(lambda_grid), digits=2), "Estimated In Red"),
           ylab="In Sample Differences")
      points(estimated_fun, type="l", col="red")
      plot(out_of_sample_dif, type="l", ylab="Out of Sample Differences",
           main=paste("MSE is ", round(MSE[i],3),  "Odd in Red"))
      lines(estimated_fun,type="l", col="red")
    }
  }
  
  if(doplot){
    par(mfrow=c(1,1))
    plot(lambda_grid, MSE)
  }
  if(plotThresh) return(Thresh_list)
  # PThresh_Chosen <- Percent_Thresholded[which.min(MSE)]
  Position <- which.min(MSE)
  opt_lambda <- lambda_grid[Position]
  Relative_Position <- Position/length(lambda_grid)
  return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda, "Position"=Relative_Position))
}

Easy_SP <- function(x, y, ..., L2=TRUE){
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  x <- sort(x)
  n <- length(x)
  i <- 1:n
  sp <- (2/pi) * asin(sqrt((i-.5)/n))
  si <- (2/pi) * asin(sqrt(y(x, ...)))
  z <- sp-si
  if(L2){
    STAT <- sum(z^2)
    return(STAT)
  }
  STAT <- sum(abs(z))
}

### only be used in package PoweR
SP_Test_Uniformity <- function(x){
  # if(min(x) < 0 || max(x) > 1) stop("Data must be b/w 0 and 1")
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  x <- sort(x)
  n <- length(x)
  i <- 1:n
  sp <- (2/pi) * asin(sqrt((i-.5)/n))
  si <- (2/pi) * asin(sqrt(punif(x)))
  z <- sp-si
  STAT <- sum(z^2)
  RVAL <- list(statistic = c(CO=STAT), p.value = 0, 
               method = "Orndorff uniformity test", data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

SP_Test<- function(x, y, ...,  Chosen_Threshold, 
                   doplot=F, includescale=F, minthreshlevel=2,
                   wf="haar", diag=F, multiple=F){
  
  
  # includescale: When doing thresholding, do we include the lowest level scaling
  # coefficient? This may not make much sense if minthreshlevel is not equal to 1
  
  # minthreshlevel: What detail level do we want to threshold to? Acceptable values
  # are from 0 to nlevelsWT()-1
  
  # wf: wavelet basis. Acceptable options for now are "haar" and "la8"
  
  # diag breaks the function for use in other functions
  # Wrapper function.
  # uses wd for convenience. Change this later to use other wavelet basis
  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  
  # Use lambda, threshold appropriately
  x <- sort(x)
  n <- length(x)
  i <- 1:n
  sp <- (2/pi) * asin(sqrt((i-.5)/n))
  si <- (2/pi) * asin(sqrt(y(x, ...)))

  if(doplot){
    par(mfrow=c(2,1))
    plot(ecdf(x))
    x1<- seq(min(x),max(x), length.out=500)
    lines(x1, y=pnorm(x1))
    par(mfrow=c(2,1))
    plot(sp,si, xlim=c(0,1), ylim=c(0,1), main="SP-Plot")
    abline(a=0, b=1)
    plot(sp-si)
  }
  z <- sp-si

  library(wavethresh)
  if(wf=="haar") model_wd <- wd(z, filter.number=1, family="DaubExPhase")
  if(wf=="la8") model_wd <- wd(z, filter.number=8, family="DaubLeAsymm")

  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=minthreshlevel:(Largest_Level-1)))

  if(includescale){
    scaling <-scaling <- accessC(model_wd, level=0)
  }else{
    scaling <- NULL
  }
  coefs <- c(details, scaling) # The vector of wavelet coefficients
  # if includescale=T, then it includes the scaling coefficient
  # Only includes detail coefficients geq to minthreshlevel
  
  thresholded_details <- sapply(coefs, function(x){
    if(abs(x) < Chosen_Threshold) x <- 0
    x
  })
  
  #### Remaining coefficients
  if(minthreshlevel >= 1){
    remaining_details <- unlist(testfun(model_wd, level=0:(minthreshlevel-1)))
  }else if(minthreshlevel==0){
    remaining_details <- NULL
  }else{
    stop("Something has gone wrong with minthreshlevel")
  }
  if(!includescale) remaining_details <- append(remaining_details, accessC(model_wd, level=0))
  # option to return coefficients
  # if 
  
  if(multiple){ 
    print("This returns all the wavelet coefficients -- when option includescale is TRUE, it returns 
the scaling coefficient as well. When includescale is FALSE, it does not")
    return(c(coefs, remaining_details))
  }
  
  threshold <- c(thresholded_details, remaining_details) 
  
  STAT <- sum(threshold^2)

  if(diag){
    thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                             value=Chosen_Threshold, levels=minthreshlevel:(nlevelsWT(model_wd)-1))

    details <- unlist(testfun(thresh_test, level=minthreshlevel:(Largest_Level-1)))
    Percent_Thresholded <- sum(details == 0) /length(details)
    reconstructed_dif <- wr(thresh_test)
    par(mfrow=c(1,2))
    plot(z, main="Actual", ylab="sp-si")
    plot(reconstructed_dif, main="Reconstructed")
    plot(z, main="Actual", type="l", ylab="sp-si")
    plot(reconstructed_dif, main="Reconstructed", type="l")

    return(list("Test Statistic"=STAT, "Percent of Eligible Coefficients Thresholded"=Percent_Thresholded)) 
  }
  STAT
}
SP_Two_Easy <- function(x, y){
  x <- sort(x)
  y <- sort(y)
  ex <- ecdf(x)
  ey <- ecdf(y)
  lenx <- length(x)
  leny <- length(y)
  sxdif  <- c(asin(sqrt((1:lenx)/lenx)) - asin(sqrt(ey(x))), asin(sqrt(ex(y))) - asin(sqrt((1:leny)/leny)))
  STAT <- sum(abs(sxdif)^2)*(lenx * leny)/(lenx + leny)^2 #L2 #Normalized
  return(STAT)
}

# ONLY FOR CONTINUOUS DISTRIBUTIONS
SP_Two_Quick <- function(x, y, L2=TRUE){
  # Delete Later
  x <- sort.int(x)
  y <- sort.int(y)
  lenx <- length(x)
  leny <- length(y)
  com <- c(x, y)
  rcom <- rank(com)
  i <- 1:lenx
  j <- 1:leny
  dif1 <- asin(sqrt(i/lenx)) - asin(sqrt((rcom[1:lenx] - i) / leny))
  dif2 <- asin(sqrt((rcom[(lenx + 1) : (lenx + leny)] - j) / lenx)) - asin(sqrt((j/leny)))
  sxdif <- c(dif1, dif2)
  if(L2){ 
    STAT <- sum(sxdif^2)*(lenx * leny)/(lenx + leny)^2 
    return(STAT)
  }
  STAT <- sum(abs(sxdif))*(lenx * leny)/(lenx + leny)^2 
  return(STAT)
}



