
Even_Odd_Grid_After <- function(x, y, ...,
                                grid_fine=10, dyadic_after=TRUE, 
                                peek_scaling=TRUE, plotEstFun=F, doplot=F,
                                debug=F, diag=F){ 
  # Use dyadic_after if the test wil transform to a dyadic sample size later
  
  # grid_fine: How fine the cross validation grid is
  # dyadic_after: Making the lambda value smaller or larger 
  # peek_scaling: Effects our scaling. Do we look at all the data when we scale,
  # or do we just look at the in sample data
  # If true, we do, if false we don't
  
  # This version peaks at the even and odd limits and expands the function appropriately
  # The previous version of the function applies the gridding with all the data
  # The previous version gridded even seperately from odd
  require(wavethresh)
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  
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
  
  # Getting the difference between the ecdf(X) and its CDF
  lenx <- length(x)
  x <- sort(x)
  F.x <- ecdf(x)
  F.x <- F.x(x)
  F.x <- F.x - (.5 /lenx)
  Dif_X <- F.x - y(x, ...) # remove
  
  # Split Samples
  
  # Split the Samples
  if(lenx %% 2==0){
    #print("its even!")
    
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    dif_even <- Dif_X[ind_even]
    data_even <- x[ind_even]
    len_even <- lenx / 2
    dif_odd <- Dif_X[ind_odd]
    data_odd <- x[ind_odd]
    len_odd <- len_even
  } else{
    #print("its odd!")
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    dif_even <- Dif_X[ind_even]
    data_even <- x[ind_even]
    len_even <- (lenx - 1) / 2
    dif_odd <- Dif_X[ind_odd]
    data_odd <- x[ind_odd]
    len_odd <- len_even + 1
  } 
  
  
  # Scaling and putting things on a grid
  # x_scale <- (x - min(x)) / (max(x)- min(x))
  
  # Here is the where the function changes
  # even_data_scaled <- (data_even - min(data_even)) / (max(data_even)- min(data_even))
  # odd_data_scaled <- (data_odd - min(data_odd)) / (max(data_odd)- min(data_odd))
  if(peek_scaling){
    even_data_scaled <- (data_even - min(x)) / (max(x)- min(x))
    odd_data_scaled <- (data_odd - min(x)) / (max(x)- min(x))
    if(doplot){
      par(mfrow=c(2,1))
      plot(even_data_scaled, dif_even, 
           main="Peaking at scaling", ylab="Fhat(x)-F(x)",
           xlim=c(0,1))
      abline( h=0)
      
      plot(odd_data_scaled, dif_odd, 
           main="Peaking at scaling", ylab="Fhat(x)-F(x)",
           xlim=c(0,1))
      abline( h=0)
      
    }
  } else{
    even_data_scaled <- (data_even - min(data_even)) / (max(data_even)- min(data_even))
    odd_data_scaled <- (data_odd - min(data_odd)) / (max(data_odd)- min(data_odd))
    if(doplot){
      par(mfrow=c(2,1))
      plot(even_data_scaled, dif_even, 
           main="No Peaking at odd data", ylab="Fhat(x)-F(x)",
           xlim=c(0,1))
      abline(h=0)
      plot(odd_data_scaled, dif_odd, 
           main="No Peaking at even data", ylab="Fhat(x)-F(x)",
           xlim=c(0,1))
      abline( h=0)
    }
  }
  
  if(debug){
    return(list("Scaled Even Data"=even_data_scaled, 
                "Scaled Odd Data"=odd_data_scaled,
                "Even Dif"=dif_even,
                "Odd Dif"= dif_odd))
  }
  

  
  # Big Assumption Here
  
  
  x_grid_e <- makegrid(t=even_data_scaled, y=dif_even)
  data_grid_e <- x_grid_e$gridt
  dif_grid_e <- x_grid_e$gridy

  
  x_grid_o <- makegrid(t=odd_data_scaled, y=dif_odd)
  data_grid_o <- x_grid_o$gridt
  dif_grid_o <- x_grid_o$gridy
  
  if(doplot){
    par(mfrow=c(2,1))
    plot(data_grid_e, dif_grid_e, 
         main="Peaking at odd data for scaling", ylab="Fhat(x)-F(x)",
         xlim=c(0,1))
    plot(data_grid_o, dif_grid_o, 
         main="Peaking at even data for scaling", ylab="Fhat(x)-F(x)",
         xlim=c(0,1))
  }
  

  
  
  Odd_OOS <- Find_Lambda(in_sample_dif = dif_grid_e, out_of_sample_dif = dif_grid_o,
                         in_sample_data = data_grid_e, out_of_sample_data = data_grid_o,
                         grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot)
  Odd_MSE <- Odd_OOS[["MSE Vector"]]
  Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
  
  Even_OOS <- Find_Lambda(in_sample_dif = dif_grid_o, out_of_sample_dif =dif_grid_e,
                          in_sample_data = data_grid_o, out_of_sample_data = data_grid_e,
                          grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot)
  Even_MSE <- Even_OOS[["MSE Vector"]]
  Even_Lambda <- Even_OOS[["Optimal Lambda"]]
  
  Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
  # This Chosen Lambda value needs to be adjusted to 
  # account for the fact that the even and odd samples
  # are of sample size n/2
  
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
                "Sample Size Adjusted Threshold"=Adjusted_Lambda))
  }
  return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
              "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
              "Chosen Threshold (Lambda)"= Chosen_Lambda,
              "Sample Size Adjusted Threshold"=Adjusted_Lambda))
  
}



Threshold_Test <- function(x, y, ..., grid=10, Chosen_Threshold, spacing=F){
  # Wrapper function.
  # uses wd for convenience. Change this later to use other wavelet basis
  
  # Carry out the thresholding, using a supplied value of 
  # Chosen_Threshold
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
  F.x <- ecdf(x)
  lenx <- length(x)
  n <- 2^floor(log2(lenx))
  z <- seq(min(x), max(x), length.out=n)
#   library(waveslim)
#   F.dwt <- dwt(F.x(z) - y(z))
#   coefs <- unlist(F.dwt)
  if(spacing){
    library(wavethresh)
    F.x <- ecdf(x)
    model_wd <- wd(F.x(x) - y(x, ...),  filter.number=8, family="DaubLeAsymm")
  }else{
    model_wd <- wd(F.x(z) - y(z, ...),  filter.number=8, family="DaubLeAsymm")
  }


  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=1:(Largest_Level-1)))
  scaling <- accessC(model_wd, level=1)
  coefs <- c(details, scaling)
  
  threshold <- sapply(coefs, function(x){
    if(abs(x) < Chosen_Threshold) x <- 0
    x
  })
  STAT <- sum(threshold^2)
  STAT
}



Find_Lambda <- function(in_sample_dif, out_of_sample_dif, 
                        in_sample_data, out_of_sample_data,
                        grid_fine, plotEstFun, doplot){
  # internal function for cross validation. Not to be used 
  
  # Input: A 'in sample' x and y and an out of sample x and y
  # Output: A list containing the MSE's used as well as the optimal lambda value found
  
  
  # Decompose input function
  model_wd <- wd(in_sample_dif, filter.number=8, family="DaubLeAsymm")
  
  # Make lambda grid
  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=1:(Largest_Level-1)))
  scaling <- accessC(model_wd, level=1)
  Largest_Coefficient <- max(abs(c(details, scaling)))
  Smallest_Coefficient <- min(abs(c(details, scaling)))
  lambda_grid <- seq(Smallest_Coefficient, 
                     Largest_Coefficient, length.out=grid_fine)
  
  # Making sure last threshold thresholds everything
  lambda_grid <- append(lambda_grid, Largest_Coefficient + .1*sd(lambda_grid))
  
  # For values in the lambda grid,
  # Perform thresholding, reconstruction, and MSE comparison
  MSE <- vector(length=length(lambda_grid))
  for(i in seq_along(lambda_grid)){
    # Thresholding
    thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                             value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1))#, verbose=TRUE)
    
    # Reconstruction
    reconstructed_fun <- wr(thresh_test)
    
    # Here, we estimate the odd from the even, and vice versa
    
    # First, we need to do some endpoint correction
    # Then, we take our estimate as the midpoint between entries
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    if(min(in_sample_data) > min(out_of_sample_data)){
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
    
    if(plotEstFun){
      #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
      par(mfrow=c(2,1))
      plot(in_sample_dif, type="l",
           main=paste("Lambda Postion=", i/length(lambda_grid)))
      points(estimated_fun, type="l", col="red")
      plot(out_of_sample_dif, type="l", main=paste("MSE is ", round(MSE[i],3),  "Odd in Red"))
      lines(estimated_fun,type="l", col="red")
    }
    
    
  }
  if(doplot){
    par(mfrow=c(1,1))
    plot(lambda_grid, MSE)
  }
  opt_lambda <- lambda_grid[which.min(MSE)]
  return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
}



# 
# Even_Odd_Threshold2 <- function(x, y, ..., lambda=FALSE, doplot=T, 
#                                grid_fine=10, opt=NULL){ 
#   # This version of the function applies the gridding with all the data
#   # The previous version gridded even seperately from odd
#   
#   if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
#   # x <- rnorm(40)
#   
#   lenx <- length(x)
#   x <- sort(x)
#   # One Sample
#   if (is.list(y)) 
#     y <- names(y)
#   if (is.function(y)) 
#     funname <- as.character(substitute(y))
#   if (is.character(y)) 
#     funname <- y
#   y <- get(funname, mode = "function", envir = parent.frame())
#   if (!is.function(y)) 
#     stop("'y' must be numeric or a function or a string naming a valid function")
# 
#   # Getting the difference between the ecdf(X) and its CDF
#   lenx <- length(x)
#   x <- sort(x)
# 
#   if(lenx %% 2==0){
#     print("its even!")
#     ind <- 1:lenx
#     ind_even <- seq(2, lenx, by=2)
#     ind_odd <- seq(1, lenx-1, by=2)
#     x_even <- x[ind_even]
#     len_even <- lenx / 2
#     x_odd <- x[ind_odd]
#     len_odd <- len_even
#   } else{
#     print("its odd!")
#     ind <- 1:lenx
#     ind_even <- seq(2, lenx, by=2)
#     ind_odd <- seq(1, lenx-1, by=2)
#     x_even <- x[ind_even]
#     len_even <- (lenx - 1) / 2
#     x_odd <- x[ind_odd]
#     len_odd <- len_even + 1
#   } 
#   lenx <- length(x)
#   x <- sort(x)
#   F.x <- ecdf(x)
#   F.x <- F.x(x)
#   F.x <- F.x - (.5 /lenx)
#   Dif_X <- F.x - y(x) # remove
#   length(Dif_X)
#   
#   # Scaling and putting things on a grid
#   x_scale <- (x - min(x)) / (max(x)- min(x))
#   
#   x_grid_e <- makegrid(t=x_scale[ind_even], y=Dif_X[ind_even])
#   data_grid_e <- x_grid_e$gridt
#   dif_grid_e <- x_grid_e$gridy
#   
#   x_grid_o <- makegrid(t=x_scale[ind_odd], y=Dif_X[ind_odd])
#   data_grid_o <- x_grid_o$gridt
#   dif_grid_o <- x_grid_o$gridy
#   
# #   # even
# #   F.x <- ecdf(x_even)
# #   F.x <- F.x(x_even)
# #   F.x <- F.x - (.5 /len_even)
# #   Dif_Even <- F.x - y(x_even, ...)
# #   #odd
# #   F.x <- ecdf(x_odd)
# #   F.x <- F.x(x_odd)
# #   F.x <- F.x - (.5 /len_odd)
# #   Dif_Odd <- F.x - y(x_odd, ...)
#   
# 
#   
#   Find_Lambda <- function(in_sample_dif, out_of_sample_dif, in_sample_data, out_of_sample_data){
#     # in_sample_dif = even_grid
#     #out_of_sample_dif= odd_grid
#     len_in_samp <- length(in_sample_dif)
#     
#     # Do wavelet decomp on even entries
#     model_wd <- wd(in_sample_dif)#, family="DaubExPhase", filter.number=1) # Add more options here later
#     
#     # Extract largest wavelet coefficient
#     Largest_Level <- nlevelsWT(model_wd)
#     testfun <- Vectorize(accessD.wd, "level")
#     details <- unlist(testfun(model_wd, level=2:(Largest_Level-1)))
#     scaling <- accessC(model_wd, level=2)
#     Largest_Coefficient <- max(abs(c(details, scaling)))
#     Smallest_Coefficient <- min(abs(c(details, scaling)))
#     lambda_grid <- seq(Smallest_Coefficient, 
#                        Largest_Coefficient, length.out=grid_fine)
#     
#     # Thresholding
#     MSE <- vector(length=length(lambda_grid))
#     for(i in seq_along(lambda_grid)){
#       #Threshold
#       thresh_test <- threshold(model_wd, type="hard", policy="manual", 
#                                value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1))#, verbose=TRUE)
#       
#       # Reconstruct
#       reconstructed_fun <- wr(thresh_test)
#       #par(mfrow=c(2,2))
#       # HERE WE DO THE INTERPOLATION
#       
#       # Need to do estimate the odd function with the even function
#       ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
#       if(min(in_sample_data) > min(out_of_sample_data)){
#         # Means in sample is 'even' and 
#         # out of sample is 'odd'
#         # Need to estimate the first data dif
#         Dif_First <- reconstructed_fun[1] / 2 # Take midpoint  
#         estimated_fun <- ma(reconstructed_fun)
#         estimated_fun <- c(Dif_First, head(estimated_fun, -1))
#         #reconstructed_fun <- c(Dif_First, ma(reconstructed_fun))
#         
#         # A potentially more complicated option
#         # library(truncdist)
#         # Interpolate minimum
#         # begin_point <- extrunc("norm", b=min(in_sample_data), ...)
#         # print(paste("Begin Point is", begin_point))
#         # print(min(in_sample_data))
#         # if_End <- ((1 - .05)/ len_even) - y(begin_point)
#         #print(Dif_End)
#       }else{
#         #Interpolate maximum
#         Dif_End <- tail(reconstructed_fun , 1) / 2
#         estimated_fun <- ma(reconstructed_fun)
#         estimated_fun <- c(head(estimated_fun, -1), Dif_End)
#         #Estimated odd end point
#         # end_point <- extrunc("norm", a=max(in_sample_data), ...)
#         # print(paste("End Point is", end_point))
#         #Dif_End <- ((len_even-.5)/ len_even) - y(end_point)
#       }
#       
#       
#       errors <- estimated_fun - out_of_sample_dif
#       
#       
#       MSE[i] <- sum(errors^2)
#       if(doplot){
#         layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
#         plot(in_sample_dif, type="l")
#         plot(estimated_fun, type="l")
#         plot(out_of_sample_dif, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
#         lines(estimated_fun,type="l", col="red")
#       }
#       
#       
#     }
#     plot(lambda_grid, MSE)
#     opt_lambda <- lambda_grid[which.min(MSE)]
#     return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
#   }
#   print(length(Dif_Even))
#   Odd_OOS <- Find_Lambda(in_sample_dif = dif_grid_e, out_of_sample_dif = dif_grid_o,
#                          in_sample_data =   data_grid_e, out_of_sample_data = data_grid_o)
#   Odd_MSE <- Odd_OOS[["MSE Vector"]]
#   Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
#   Even_OOS <- Find_Lambda(in_sample_dif = dif_grid_o, out_of_sample_dif =dif_grid_e,
#                           in_sample_data = data_grid_o, out_of_sample_data = data_grid_e)
#   Even_MSE <- Even_OOS[["MSE Vector"]]
#   Even_Lambda <- Even_OOS[["Optimal Lambda"]]
#   Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
#   print(Chosen_Lambda)
#   return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
#               "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
#               "Chosen Threshold (Lambda)"= Chosen_Lambda))
#   
#   
#   #   
#   #   # Do wavelet decomp on even entries
#   #   even_wd <- wd(even_grid)#, family="DaubExPhase", filter.number=1) # Add more options here later
#   #   
#   #   # Find Universal Threshold
#   #   wavelet_coefs <- NULL
#   #   for(i in 3:(nlevelsWT(even_wd)-1))  
#   #     wavelet_coefs <- c(wavelet_coefs, accessD(even_wd, level=i))
#   #   noise_level <- mad(wavelet_coefs) #noise.level
#   #   universal_threshold <- noise_level * sqrt(2*log(even_len))
#   #   
#   #     # Create lambda grid
#   #   min_lambda <- universal_threshold - 15 * sd(even_grid)
#   #   max_lambda <- universal_threshold + 3 * sd(even_grid)
#   #   print(min_lambda)
#   #   print(max_lambda)
#   #   if(min_lambda < 0) min_lambda=0
#   #   lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
#   #   
#   #   MSE <- vector(length=length(lambda_grid))
#   #   for(i in seq_along(lambda_grid)){
#   #     #Threshold
#   #     thresh_test <- threshold(even_wd, type="hard", policy="manual", 
#   #                              value=lambda_grid[i], levels=0:(nlevelsWT(even_wd)-1),
#   #                              verbose=TRUE)
#   #     #print(lambda_grid[i])
#   #     #print(thresh_test)
#   #     # Reconstruct
#   #     reconstructed_fun <- wr(thresh_test)
#   #     #par(mfrow=c(2,2))
#   #     errors <- reconstructed_fun - odd_grid
#   #     
#   #     MSE[i] <- sum(errors^2)
#   #     layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
#   #     plot(even_grid, type="l")
#   #     plot(reconstructed_fun, type="l")
#   #     plot(odd_grid, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
#   #     lines(reconstructed_fun,type="l", col="red")
#   #   }
#   #   
#   #   if(doplot){
#   #     par(mfrow=c(1,2))
#   #     plot(x_odd, Dif_Odd, type="l")
#   #     plot(x_even, Dif_Even, type="l")
#   #   }
#   #   
#   #   return(MSE)
#   
#   # return(list("Even X"=x_even, "Odd X"=x_odd, "Even Dif"=Dif_Even, "Odd Dif"=Dif_Odd))
# }
# 
# 
# 
# ### OLD CODE
# 
# 
# 
# 
# 
# 
# # # Transforming to grid (temporary solution)
# # x_even_scale <- (x_even - min(x_even)) / (max(x_even) - min(x_even))
# # x_odd_scale <- (x_odd - min(x_odd)) / (max(x_odd) - min(x_odd))
# # 
# # even_grid <- makegrid(t=x_even_scale, y=Dif_Even)
# # even_grid <- even_grid$gridy
# # 
# # odd_grid <- makegrid(t=x_odd_scale, y=Dif_Odd)
# # odd_grid <- odd_grid$gridy
# # 
# # Find_Lambda <- function(in_sample_data, out_of_sample_dif){
# #   # in_sample_data = even_grid
# #   #out_of_sample_dif= odd_grid
# #   
# #   # Do wavelet decomp on even entries
# #   model_wd <- wd(in_sample_data)#, family="DaubExPhase", filter.number=1) # Add more options here later
# #   
# #   # model_wd = even_wd
# #   # Find Universal Threshold
# #   wavelet_coefs <- NULL
# #   for(i in 3:(nlevelsWT(model_wd)-1))  
# #     wavelet_coefs <- c(wavelet_coefs, accessD(model_wd, level=i))
# #   noise_level <- mad(wavelet_coefs) #noise.level
# #   universal_threshold <- noise_level * sqrt(2*log(even_len))
# #   
# #   # Create lambda grid
# #   min_lambda <- universal_threshold - 15 * sd(in_sample_data)
# #   max_lambda <- universal_threshold + 3 * sd(in_sample_data)
# #   if(min_lambda < 0) min_lambda=0
# #   lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
# #   
# #   # Thresholding
# #   MSE <- vector(length=length(lambda_grid))
# #   for(i in seq_along(lambda_grid)){
# #     #Threshold
# #     thresh_test <- threshold(model_wd, type="hard", policy="manual", 
# #                              value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1),
# #                              verbose=TRUE)
# #     
# #     # Reconstruct
# #     reconstructed_fun <- wr(thresh_test)
# #     #par(mfrow=c(2,2))
# #     errors <- reconstructed_fun - out_of_sample_dif
# #     
# #     MSE[i] <- sum(errors^2)
# #     if(doplot){
# #       layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
# #       plot(even_grid, type="l")
# #       plot(reconstructed_fun, type="l")
# #       plot(odd_grid, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
# #       lines(reconstructed_fun,type="l", col="red")
# #     }
# #   }
# #   opt_lambda <- lambda_grid[which.min(MSE)]
# #   return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
# # }
# # 
# # Odd_OOS <- Find_Lambda(in_sample_data = even_grid, out_of_sample_dif = odd_grid)
# # Odd_MSE <- Odd_OOS[["MSE Vector"]]
# # Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
# # Even_OOS <- Find_Lambda(in_sample_data = odd_grid, out_of_sample_dif = even_grid)
# # Even_MSE <- Even_OOS[["MSE Vector"]]
# # Even_Lambda <- Even_OOS[["Optimal Lambda"]]
# # Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
