# EWMA functions (qiu is getting crowded)

# Generic function for finding control limits 
# in an EWMA framework


### The problem with this approach is that it is too slow
### We can speed it up by pre-calculating the ecdf for the in control
### data

# Solution:
# Use Switch Statements to Branch 

EWMA_Find_CL_Slow <- function(ic_data = rnorm(500),
                         lambda = .05,
                         control_limit = .12,
                         m = 5,
                         ICdist = "rnorm",
                         IC_dist_ops = NULL,
                         bootstrap_samples = 3000,
                         test_stat = "cvm") {
  # Calculate d0 via bootstrapping
  D_n <- numeric(bootstrap_samples)
  for (i in 1:bootstrap_samples) {
    sampled_data <- sample(ic_data, replace = T, size = m)
    D_n[i] <- test_stat(ic_data, sampled_data)
  }
  d0 = mean(D_n)
  # Initializing Variables
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    
    
    D_n <- test_stat(ic_data, data)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    u <- append(u, u_nK)
  }
  return(list("control_stats" = u, "Time OOC" = i))
}

EWMA_Find_CL_CVM <- function(ic_data = rnorm(500),
                              lambda = .05,
                              control_limit = .06,
                              m = 5,
                              ICdist = "rnorm",
                              IC_dist_ops = NULL,
                              bootstrap_samples = 3000,
                              track_candidates=NULL) {
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)

  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  i <- 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)
  ind <- 1:len_ic
  
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    D_n <- Fast_TS_CVM(sorted_ic_data=sorted_ic_data, new_data=data, 
                lenx=len_ic, 
                j=j, i=ind, leny=m)
    # print(D_n)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
  }
  
  # If we have a large upper bound, we can calculate all our upper bounds below that number
  if(!is.null(track_candidates)){
    track_ucl <-numeric(length(track_candidates))
    
    for(j in seq_along(track_candidates)){
      track_ucl[j] <- as.numeric(which.min(track_candidates[j] > u)) - 1 #investigate later
    }
    names(track_ucl) <- as.character(track_candidates)
    return(list("Time OOC"=i, "Lower CLs"=track_ucl))
    
  }
  return(list("control_stats" = u, "Time OOC" = i))
}

Fast_Bootstrap_CVM <- function(ic_data, m, bootstrap_samples){
  # This is inaccurate (5/5/2016)
  # Possibly because of ties
  D_n <- numeric(bootstrap_samples)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  i <- 1:lenx
  leny <- m
  j <- 1:leny
  sd_data <- sd(ic_data)
  for (k in 1:bootstrap_samples) {
    sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
    #sampled_data <- sampled_data + rnorm(m, 0, .01*sd_data)
    #sampled_data <- rnorm(5)
    #print(sampled_data)
    ranks <- rank(c(sort_ic, sampled_data))
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    D_n[k] <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
  }
  d0 = mean(D_n)
}

In_Control_Dist_CvM <- function(ic_data, m, bootstrap_samples){
  # This is inaccurate (5/5/2016)
  # Possibly because of ties
  D_n <- numeric(bootstrap_samples)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  i <- 1:lenx
  leny <- m
  j <- 1:leny
  sd_data <- sd(ic_data)
  for (k in 1:bootstrap_samples) {
    sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
    #sampled_data <- sampled_data + rnorm(m, 0, .01*sd_data)
    #sampled_data <- rnorm(5)
    #print(sampled_data)
    ranks <- rank(c(sort_ic, sampled_data))
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    D_n[k] <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
  }
  D_n
}

In_Control_Dist_CvM_Exp <- function(ic_data, m, bootstrap_samples, ICdist, ic_params){
  # This is inaccurate (5/5/2016)
  # Possibly because of ties
  D_n <- numeric(bootstrap_samples)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  i <- 1:lenx
  leny <- m
  j <- 1:leny
  sd_data <- sd(ic_data)
  
  # Newly simulated data
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  # new_data <- (lenx) # breaks upon usage of anything other than default params
  # new_data <- do.call(ic_dist, c(list(lenx), ic_params)) 
  
  
  
  for (k in 1:bootstrap_samples) {
    sampled_data <- do.call(rIC, c(list(m), ic_params))  #NOW A MISNOMER
    # sampled_data <- sort.int(sample(new_data, replace = T, size = m))
    #sampled_data <- sampled_data + rnorm(m, 0, .01*sd_data)
    #sampled_data <- rnorm(5)
    #print(sampled_data)
    ranks <- rank(c(sort_ic, sampled_data))
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    D_n[k] <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
  }
  D_n
}


# 
# Fast_Bootstrap_CVM_Ties <- function(ic_data, m, bootstrap_samples){
#   # Trying to make this compatiable with ties
#   D_n <- numeric(bootstrap_samples)
#   sort_ic <- sort(ic_data)
#   lenx <- length(ic_data)
#   rankx <- rank(sort_ic)
#   for (i in 1:bootstrap_samples) {
#     sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
#     ranks <- rank(c(sort_ic, sampled_data))
#     U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
#     D_n[i] <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
#   }
#   d0 = mean(D_n)
# }


Fast_TS_CVM <- function(sorted_ic_data, new_data, lenx, i, j, leny){
 # Internal use ONLY 
  ranks <- rank(c(sorted_ic_data, new_data))
  U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
  STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
}


EWMA_CVM_OOC <- function(ic_data=rnorm(500), lambda=.05, control_limit, m=5, exact=FALSE, 
                          tau=0, 
                          ICdist="rnorm", IC_dist_ops=NULL,
                          OOCdist="rnorm", OOC_dist_ops=NULL,
                          bootstrap_samples=3000){
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())

  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)
  ind <- 1:len_ic
  
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 
      
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    # Calculate new value of test statisic
    D_n <- Fast_TS_CVM(sorted_ic_data=sorted_ic_data, new_data=data, 
                       lenx=len_ic, 
                       j=j, i=ind, leny=m)
    # print(D_n)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
  }
  return(list("control_stats" = u, "Time OOC" = i))
}


Iterative_CL_Finding <- function(allow_param, cl_seq, f_in_control, margin_error){
  # allow_param = Allowance Parameter value
  # cl_seq = control limit sequence of values to be evaluated
  # f_in_control = What function to use
  # INTERNAL USE ONLY
  # FOR NOW, SET N2, N2_MIN, ETC. OUTSIDE OF THE FUNCTION
  sd_tol_init <- sd_tol
  time_ic <- numeric(length(cl_seq))
  Iterations <- 0
  arl_margin_error <- margin_error + 1 
  while(arl_margin_error > margin_error){
    for(j in seq_along(cl_seq)){
      N2 <- 0 # At least 300 runs
      sd_arl <- sd_tol + 2
      arl_track <- NULL
      while(sd_arl > sd_tol){
        current_run_length <- f_in_control(control_limit=cl_seq[j], 
                                           lambda=allow_param, m=5,
                                           ic_data=rdist(M))[["Time OOC"]]
        arl_track <- append(arl_track, current_run_length)
        sd_arl <- sd(arl_track)* 1.96 / sqrt(N2)
        arl_est <- mean(arl_track)
        if(is.na(sd_arl)) sd_arl <- sd_tol + 1
        N2 <- N2 + 1
        if(N2 < N2_min) sd_arl <- sd_tol + 1
        if(N2 %% 40==0){
          print(N2)
          print(arl_est)
          print(sd_arl)
        }
        time_ic[j] <- mean(arl_track)
      }
    }
    xi <- cl_seq
    yi <- time_ic
    
    m1 <- lm(xi ~ yi + I(yi^2))
    cl_pred <- predict(m1, newdata=data.frame(yi=ARL))
    
    cl_int <- predict(m1, newdata=data.frame(yi=ARL), interval="confidence")
    
    if(cl_pred < 0){
      cl_pred <- 0
      print("Your estimate was way off")
    }
    if(cl_int[2] < 0) cl_int[2] <- 0
    
    print(paste("X-Values are", xi))
    print(paste("Y-Values are", yi))
    
    print(paste("Confidence Interval for control limit is", cl_int[2], ",", cl_int[3])) 
    # Go the other way to establish confidence limits for the ARL
    m2 <- lm(yi ~ xi +I(xi^2))
    arl_pred <- predict(m2, newdata=data.frame(xi=cl_pred), interval="confidence")
    print(paste("Confidence Interval for ARL", ARL, "is", arl_pred[2],",", arl_pred[3]))
    
    # Determine margin error for average run length
    arl_lower_bound <- arl_pred[2]
    arl_upper_bound <- arl_pred[3]
    arl_margin_error <- arl_upper_bound - arl_lower_bound
    
    # Re-adjust parameters for control limit
    lcl_update <- cl_int[2]
    ucl_update <- cl_int[3]
    cl_seq_update <- seq(lcl_update, ucl_update, length.out=5) #5 is arbitrary
    # replace cl_seq with cl_seq_update
    cl_seq <- cl_seq_update
    
    Iterations <- Iterations + 1
    # Demand more specificity with each iteration
    # This is a little controversial
    sd_tol <- sd_tol_init/(Iterations + 1)
    print(paste("This is Iteration", Iterations))
    print(paste("Tolerance for ARL Calculations is now", sd_tol))
  }
  
  return(list("Predicted CL"=cl_pred, "CL Interval"=cl_int, "ARL Interval"=arl_pred,
              "Number of Iterations Required for Convergence"=Iterations))
}


CL_Finding_With_UB <- function(allow_param, cl_seq, f_in_control, sd_tol_init, Cores=4){
  # allow_param = Allowance Parameter value
  # cl_seq = control limit sequence of values to be evaluated
  # f_in_control = What function to use
  # INTERNAL USE ONLY
  # FOR NOW, SET N2, N2_MIN, ETC. OUTSIDE OF THE FUNCTION
  
  N2 <- 10000
  N2_min <- 500
  
  
  upper.bound.ucl <- max(cl_seq)
  num_candidates <- length(cl_seq)
    Iterations <- 0
  sd_arl <- sd_tol_init + 1 
  Num_Iter <- 5 * Cores
  require(doMC)
  require(foreach)
  N2 <- 0 # At least 300 runs
  arl_track <- NULL
  while(sd_arl > sd_tol_init){
    print(paste("sd_arl was", sd_arl, "sd_tol_init was", sd_tol_init))
    par_obj <- foreach(n=1:Num_Iter) %dopar%{
      RLs <- f_in_control(control_limit=max(cl_seq), 
                          lambda=allow_param, m=5,
                          ic_data=rdist(M),
                          track_candidates = cl_seq)
      current_run_length <- RLs[["Time OOC"]]
      names(current_run_length) <- as.character(upper.bound.ucl)
      lower_cls <- RLs[["Lower CLs"]]
      all_RLs <- c(current_run_length, lower_cls)
    }
    par_obj = matrix(unlist(par_obj), nrow=num_candidates+1, byrow=FALSE)
    arl_track <- cbind(arl_track, par_obj) # arl_track is now a matrix
    N2 <- N2 + Num_Iter
    sd_arl <- sd(arl_track[1, ])* 1.96 / sqrt(N2) # First row contains most relevent info
    arl_est <- mean(arl_track[1, ] ) 
    
    # Old Vector Code
    # par_obj <- unlist(par_obj)
    # print(par_obj)
    # arl_track <- append(arl_track, par_obj) # Keeps track of our OOC times at the current UCL
    # N2 <- N2 + Num_Iter
    # sd_arl <- sd(arl_track)* 1.96 / sqrt(N2)
    # arl_est <- mean(arl_track) # Takes the mean of the OOC times
    
    if(is.na(sd_arl)) sd_arl <- sd_tol_init + 1
    
    if(N2 < N2_min) sd_arl <- sd_tol_init + 1
    if(N2 > N2_max) sd_arl <- sd_tol_init - 1 # goes to far, just need to cut it off
    if(N2 %% 40==0){
      print(N2)
      print(arl_est)
      print(sd_arl)
      #print(arl_track)
      #time_ic[j] <- mean(arl_track)
      
    }
  }
  rownames(arl_track) <- as.character(c(upper.bound.ucl,cl_seq))
  return(arl_track)
}

Iterative_CL_Finding_MC <- function(allow_param, cl_seq, f_in_control, margin_error, Cores=4){
  # allow_param = Allowance Parameter value
  # cl_seq = control limit sequence of values to be evaluated
  # f_in_control = What function to use
  # INTERNAL USE ONLY
  # FOR NOW, SET N2, N2_MIN, ETC. OUTSIDE OF THE FUNCTION
  sd_tol_init <- sd_tol
  time_ic <- numeric(length(cl_seq))
  Iterations <- 0
  arl_margin_error <- margin_error + 1 
  Num_Iter <- 5 * Cores
  require(doMC)
  require(foreach)
  
  while(arl_margin_error > margin_error){
    for(j in seq_along(cl_seq)){
      N2 <- 0 # At least 300 runs
      sd_arl <- sd_tol + 2
      arl_track <- NULL
      while(sd_arl > sd_tol){
        par_obj <- foreach(n=1:Num_Iter) %dopar%{
          current_run_length <- f_in_control(control_limit=cl_seq[j], 
                                           lambda=allow_param, m=5,
                                           ic_data=rdist(M))[["Time OOC"]]
        }
        par_obj <- unlist(par_obj)
        print(par_obj)
        arl_track <- append(arl_track, par_obj)
        N2 <- N2 + Num_Iter
        sd_arl <- sd(arl_track)* 1.96 / sqrt(N2)
        arl_est <- mean(arl_track)
        if(is.na(sd_arl)) sd_arl <- sd_tol + 1

        if(N2 < N2_min) sd_arl <- sd_tol + 1
        if(N2 > N2_max) sd_arl <- sd_tol - 1 # goes to far, just need to cut it off
        if(N2 %% 40==0){
          print(N2)
          print(arl_est)
          print(sd_arl)
        }
        time_ic[j] <- mean(arl_track)
      }
    }
    xi <- cl_seq
    yi <- time_ic
    
    m1 <- lm(xi ~ yi + I(yi^2))
    cl_pred <- predict(m1, newdata=data.frame(yi=ARL))
    
    cl_int <- predict(m1, newdata=data.frame(yi=ARL), interval="confidence")
    
    if(cl_pred < 0){
      cl_pred <- 0
      print("Your estimate was way off")
    }
    if(cl_int[2] < 0) cl_int[2] <- 0
    
    print(paste("X-Values are", xi))
    print(paste("Y-Values are", yi))
    
    print(paste("Confidence Interval for control limit is", cl_int[2], ",", cl_int[3])) 
    # Go the other way to establish confidence limits for the ARL
    m2 <- lm(yi ~ xi +I(xi^2))
    arl_pred <- predict(m2, newdata=data.frame(xi=cl_pred), interval="confidence")
    print(paste("Confidence Interval for ARL", ARL, "is", arl_pred[2],",", arl_pred[3]))
    
    # Determine margin error for average run length
    arl_lower_bound <- arl_pred[2]
    arl_upper_bound <- arl_pred[3]
    arl_margin_error <- arl_upper_bound - arl_lower_bound
    
    # Re-adjust parameters for control limit
    lcl_update <- cl_int[2]
    ucl_update <- cl_int[3]
    cl_seq_update <- seq(lcl_update, ucl_update, length.out=5) #5 is arbitrary
    # replace cl_seq with cl_seq_update
    cl_seq <- cl_seq_update
    
    Iterations <- Iterations + 1
    # Demand more specificity with each iteration
    # This is a little controversial
    print(paste("This is Iteration", Iterations))
    sd_tol <- sd_tol_init/(Iterations + 1)
    
    if(arl_margin_error < margin_error){
      print(paste("Done! Margin of Error was", arl_margin_error, "which was less than", margin_error))
    }else{
    print(paste("Not Done yet! MoE was", arl_margin_error, "which was greater than", margin_error))
    print(paste("Tolerance for ARL Calculations is now", sd_tol))
    }
  }
  
  return(list("Predicted CL"=cl_pred, "CL Interval"=cl_int, "ARL Interval"=arl_pred,
              "Number of Iterations Required for Convergence"=Iterations))
}

Initial_CL_Approx_CVM <- function(num_candidates = 30,
                                  ARL = 200,
                                  bootstrap_samples = 5000,
                                  ic_sample_size = 2000,
                                  allowance_parameter = .05,
                                  ic_dist="rnorm",
                                  ic_params=list(mean=0, sd=1),
                                  runs = 1000,
                                  avg_runs = 10,
                                  est_high=TRUE,
                                  over_estimation_param=1.3,
                                  m=5){
  
  
  # TO DO: 
  #  CHANGE UPPER RUN LIMIT DEPENDING ON ARL
  num_allowance_parameters <- length(allowance_parameter)
  if(est_high){
    cl_est <- vector(mode="list", length=2)
    cl_est[[1]] <- matrix(nrow= avg_runs, ncol=num_allowance_parameters)
    cl_est[[2]] <- matrix(nrow= avg_runs, ncol=num_allowance_parameters)
  }else{
    cl_est <- matrix(nrow= avg_runs, ncol=num_allowance_parameters)
  }
  # cl_est <- numeric(avg_runs)
  ic_dist <- get(ic_dist, mode = "function", envir = parent.frame())
  
  for(count_i in 1:avg_runs){
    ic_data <- do.call(ic_dist, c(list(ic_sample_size), ic_params)) 
    
    # bootstrap sample from that in control data to get distribution of test statistic
    in_control_tstat <- In_Control_Dist_CvM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
    # get mean of the distribution
    mean_ic_tstat <- mean(in_control_tstat)
    
    # Upper run limit depends on what ARL you want to estimate
    p <- (1/ARL)
    sd_p <- sqrt((1-p)/p^2)
    upper_run_limit <- floor(ARL + 15 * sd_p) # yields 4000 when ARL=200
    # 10000 when ARL=500 (Could use some testing later)
    
    # number of runs is ad hoc right now
    
    # Pre-allocating
    # x <- numeric(upper_run_limit)
    
    
    
    matRuns <- matrix(nrow=runs, ncol=upper_run_limit)
    liRuns <- vector(mode="list", length=runs)
    # Get control limit estimates by using the bootstrapped in control distribution
    # for the test statistic
    
    for(i in 1:runs){
      x <- matrix(data= 0,nrow=num_allowance_parameters, ncol=upper_run_limit)
      tstat_chosen <- sample(in_control_tstat, upper_run_limit) # Faster to do this earlier
      for(j in 2:upper_run_limit){
        #tstat_chosen <- sample(in_control_tstat, 1)
        x[, j] <- allowance_parameter * (tstat_chosen[j] - mean_ic_tstat) + (1 - allowance_parameter) * x[, j - 1]
      }
      liRuns[[i]] <- x
    }
    
    maxes <- lapply(liRuns, function(x) apply(x, 1, max))
    # convert to matrix nrow= # num allowance parameters and ncol= # of runs
    matRunsLambda <- do.call(rbind, maxes)
    min_maxes <- apply(matRunsLambda, 2, min)
    
    # maxes <- apply(matRuns, 1, max)
    #min_max <- min(maxes) # This is our jumpoff point for estimating control limits
    lcl <- rep(0, num_allowance_parameters)
    ucl <- min_maxes
    # Later, we will need a lot of error checking here
    #cl_cand <- seq(lcl, ucl, length.out=num_candidates)
    
    # This will break unless lcl = 0
    cl_cand <- sapply(ucl, function(x) seq(0, x, length.out=num_candidates))
    # Generates two sequences 
    
    arl_candidates <- vector(mode="list", length=num_allowance_parameters)
    for(i in 1:num_allowance_parameters){
      cc <- matrix(nrow = runs, ncol=num_candidates)
      for(j in 1:runs){
        cc[j, ] <- sapply(cl_cand[, i] , function(cl_cand) min(which(liRuns[[j]][i, ] >= cl_cand)))
      }
      arl_candidates[[i]] <- colMeans(cc)
    }
    
    if(est_high){
      for(i in 1:num_allowance_parameters){
        f_in <- approxfun(arl_candidates[[i]], cl_cand[, i])
        f <- Vectorize(f_in)
        control_limit_estimate <- f(c(ARL, over_estimation_param*ARL))
        for(j in 1:2){
          cl_est[[j]][count_i, i] <- control_limit_estimate[j]
          #cl_est[count_i, i] <- control_limit_estimate[j]
        }
      }
    }else{
      for(i in 1:num_allowance_parameters){
        f_in <- approxfun(arl_candidates[[i]], cl_cand[, i])
        control_limit_estimate <- f_in(ARL)
        cl_est[count_i, i] <- control_limit_estimate
      }
    }
    print(count_i/avg_runs)
  }
  
  if(est_high){
    cl_est_mean <- lapply(cl_est, colMeans)
    cl_est_mean <- do.call(rbind, cl_est_mean)
    rownames(cl_est_mean) <- c(ARL, over_estimation_param*ARL)
    colnames(cl_est_mean) <- allowance_parameter
    return(cl_est_mean)
  }else{
    cl_est_mean <- colMeans(cl_est)
    names(cl_est_mean) <- allowance_parameter
    return(cl_est_mean)
  }
  
}

# 
# # Stored 5/25/2016 
# newcl_cvm <- function(num_candidates=30,
#                       ARL = 200,
#                       bootstrap_samples = 5000,
#                       ic_sample_size = 2000,
#                       allowance_parameter = .05,
#                       ic_dist="rnorm",
#                       ic_params=list(mean=0, sd=1),
#                       runs = 1000,
#                       avg_runs = 14){
#   # TO DO: 
#   #  CHANGE UPPER RUN LIMIT DEPENDING ON ARL
#   cl_est <- numeric(avg_runs)
#   ic_dist <- get(ic_dist, mode = "function", envir = parent.frame())
#   
#   for(count_i in 1:avg_runs){
#     ic_data <- do.call(ic_dist, c(list(ic_sample_size), ic_params)) 
#     
#     # bootstrap sample from that in control data to get distribution of test statistic
#     in_control_tstat <- In_Control_Dist_CvM(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples)
#     # get mean of the distribution
#     mean_ic_tstat <- mean(in_control_tstat)
#     
#     # Upper run limit depends on what ARL you want to estimate
#     p <- (1/ARL)
#     sd_p <- sqrt((1-p)/p^2)
#     upper_run_limit <- ARL + 15 * sd_p # yields 4000 when ARL=200
#     # 10000 when ARL=500 (Could use some testing later)
#     
#     # number of runs is ad hoc right now
#     
#     # Pre-allocating
#     x <- numeric(upper_run_limit)
#     cc <- matrix(nrow = runs, ncol=num_candidates)
#     matRuns <- matrix(nrow=runs, ncol=upper_run_limit)
#     
#     # Get control limit estimates by using the bootstrapped in control distribution
#     # for the test statistic
#     
#     for(i in 1:runs){
#       x <- numeric(upper_run_limit)
#       tstat_chosen <- sample(in_control_tstat, upper_run_limit) # Faster to do this earlier
#       for(j in 2:upper_run_limit){
#         #tstat_chosen <- sample(in_control_tstat, 1)
#         x[j] <- allowance_parameter * (tstat_chosen[j] - mean_ic_tstat) + (1 - allowance_parameter)*x[j - 1]
#       }
#       matRuns[i, ] <- x
#     }
#     maxes <- apply(matRuns, 1, max)
#     min_max <- min(maxes) # This is our jumpoff point for estimating control limits
#     lcl <- 0
#     ucl <- min_max
#     # Later, we will need a lot of error checking here
#     cl_cand <- seq(lcl, ucl, length.out=num_candidates)
#     
#     for(i in 1:runs){
#       cc[i,] <- sapply(cl_cand, function(cl_cand) min(which(matRuns[i,] >= cl_cand)))
#     }
#     
#     # averaging over all the runs
#     arl_candidates <- colMeans(cc)
#     
#     # plot(cl_cand, colMeans(cc))
#     # f <- approxfun(cl_cand,arl_candidates)
#     f_in <- approxfun(arl_candidates, cl_cand)
#     control_limit_estimate <- f_in(ARL)
#     cl_est[count_i] <- control_limit_estimate
#     print(count_i/avg_runs)
#   }
#   cl_est_mean <- mean(cl_est)
# }
# Depricated 
# Estimate_Control_Limit_EWMA_Gamma <-
#   function(M,
#            rdist=rgamma,
#            m,
#            ARL,
#            tstat,
#            allowance_parameter,
#            reps = 1000,
#            lowerARL = NULL,
#            upperARL = NULL) {
#     # Code is obviously verbose, but whatever
#     # Only one allowance parameter value allowed
#     #rdist should be a character
#     if(is.null(lowerARL)) lowerARL <- ARL - (ARL/5)
#     if(is.null(upperARL)) upperARL <- ARL + (ARL/5)
#     
#     tstat_parameters <- Estimate_Allowance_Gamma(
#       ic_dist = "rnorm",
#       M = M,
#       bootstrap_samples = 2500,
#       tstat = tstat,
#       m = m
#     )
#     
#     tshape <- tstat_parameters[1]
#     trate <- tstat_parameters[2]
#     tmean <- tstat_parameters[3]
#     if (is.character(rdist)) {
#       rdist <- get(rdist, mode = "function", envir = parent.frame())
#     }
#     #len <- ARL
#     len <- c(ARL, lowerARL, upperARL)
#     #names(len) <- c("ARL", "lowerARL", "upperARL")
#     print(len)
#     ucl_guess <- numeric(length(len))
#     names(ucl_guess) <- c("ARL", "lowerARL", "upperARL")
#     
#     for(i in seq_along(len)){
#       max_x <- numeric(reps)
#       
#      for(j in 1:reps){
#        temp_len <- len[i]
#        x <- numeric(temp_len)
#        
#        for(k in 2:temp_len){
#          temp <- rdist(1, shape=tshape, rate=trate)
# 
# 
#          x[k] <- allowance_parameter * (temp-tmean) + (1-allowance_parameter)*x[k-1]
#          #x[k] <- allowance_parameter * ( rdist(1, shape=tshape, rate=trate) - tmean)
#         #                          + (1 - allowance_parameter) * x[k-1]
#          
#        }
#        max_x[j] <- max(x)
#       }
#       ucl_guess[i] <- mean(max_x)
#     }
#     
#     return(list("Estimated Control Limit Mean"=ucl_guess["ARL"],
#                 "Estimated Control Limit UB"=ucl_guess["lowerARL"],
#                 "Estimated Control Limit LB"=ucl_guess["upperARL"]))
#     
#   }

# Depricated 
# newccfun <- function(num_candidates=30, 
#                      lcl=.03, 
#                      ucl=.055, 
#                      bootstrap_samples=5000, 
#                      ic_sample_size=2000,
#                      allowance_parameter=.05){
#   # Creating in control data
#   ic_data <- rnorm(ic_sample_size)
#   # bootstrap sample from that in control data to get distribution of test statistic
#   in_control_tstat <- In_Control_Dist_CvM(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples)
#   # get mean of the distribution
#   mean_ic_tstat <- mean(in_control_tstat)
#   
#   # Upper run limit depends on what ARL you want to estimate
#   upper_run_limit <- 4000 
#   runs <- 1000
#   control_limit_candidates <- seq(lcl, ucl, length.out=num_candidates)
#   
#   # Pre-allocating
#   x <- numeric(upper_run_limit)
#   cc <- matrix(nrow = runs, ncol=num_candidates)
#   
#   # Get control limit estimates by using the bootstrapped in control distribution
#   # for the test statistic
#   
#   for(i in 1:runs){
#     x <- numeric(upper_run_limit)
#     for(j in 2:upper_run_limit){
#       tstat_chosen <- sample(in_control_tstat, 1)
#       x[j] <- allowance_parameter * (tstat_chosen - mean_ic_tstat) + (1 - allowance_parameter)*x[j - 1]
#     }
#     cc[i,] <- sapply(control_limit_candidates, function(control_limit_candidates) min(which(x > control_limit_candidates)))
#     #print(i/runs)
#   }
#   
#   cc_means <- colMeans(cc)
#   names(cc_means) <- control_limit_candidates
#   
#   plot(control_limit_candidates, cc_means)
#   print("done!")
#   cc_means
# }



