# EWMA functions (qiu is getting crowded)

# Generic function for finding control limits 
# in an EWMA framework


### The problem with this approach is that it is too slow
### We can speed it up by pre-calculating the ecdf for the in control
### data

# Solution:
# Use Switch Statements to Branch 

CL_Finding_With_UB <- function(allow_param, 
                               cl_seq, 
                               f_in_control, 
                               sd_tol_init, 
                               Cores=4, 
                               max_iter=5000,
                               pval_approach=FALSE, ...){
  # allow_param = Allowance Parameter value
  # cl_seq = control limit sequence of values to be evaluated
  # f_in_control = What function to use
  # pval_approach: if using apvalue approach, smaller limits
  # FOR NOW, SET N2, N2_MIN, ETC. OUTSIDE OF THE FUNCTION
  
  N2_max <- max_iter # Run a maximum of this many times
  N2_min <- 500 # Run for at least 500 times
  
  upper.bound.ucl <- max(cl_seq)
  num_candidates <- length(cl_seq)
  
  # Parallel processing
  require(doMC)
  require(foreach)
  Num_Iter <- 5 * Cores
  
  sd_arl <- sd_tol_init + 1  # setting while loop condition

  arl_track <- NULL
  
  # what limit will results in the largest arl?
  if(pval_approach){
    choose_cl <- min(cl_seq)
  }else{
    choose_cl <- max(cl_seq)
  }
  N2 <- 0

  while(sd_arl > sd_tol_init){
    
    par_obj <- foreach(n=1:Num_Iter) %dopar%{
      RLs <- f_in_control(control_limit=choose_cl, 
                          lambda=allow_param, 
                          m=m,
                          ic_data= do.call(rdist, c(list(M), ic_params)),
                          track_candidates = cl_seq, 
                          ...)
      
      current_run_length <- RLs[["Time OOC"]]
      
      names(current_run_length) <- as.character(choose_cl)
      
      lower_cls <- RLs[["Lower CLs"]]

      if(pval_approach){
        # vector of run lengths for each control limit in sequence
        all_RLs <- c(lower_cls) # on 10/20 confirmed this is correct
      }else{
        all_RLs <- c(current_run_length, lower_cls) # Is this an error? check later
      }
    }
    if(pval_approach){
      par_obj = matrix(unlist(par_obj), nrow=num_candidates, byrow=FALSE)
      
    }else{
      par_obj = matrix(unlist(par_obj), nrow=num_candidates + 1, byrow=FALSE)
      
    }
    arl_track <- cbind(arl_track, par_obj) # arl_track is now a matrix
    N2 <- N2 + Num_Iter
    
    if(pval_approach){
      arl_est <- mean(arl_track[dim(arl_track)[1], ] )    
      sd_arl <- sd(arl_track[dim(arl_track)[1], ])* 1.96 / sqrt(N2) # First row contains most relevent info
      
    }else{
      arl_est <- mean(arl_track[1, ] )      
      sd_arl <- sd(arl_track[1, ])* 1.96 / sqrt(N2) # First row contains most relevent info
      
    }

    
    if(is.na(sd_arl)) sd_arl <- sd_tol_init + 1
    
    if(N2 < N2_min) sd_arl <- sd_tol_init + 1
    if(N2 > N2_max) sd_arl <- sd_tol_init - 1 # goes to far, just need to cut it off
    if(N2 %% 2==0){
      # temporary
      print(N2)
      print(arl_est)
      print(sd_arl)
    }
  }
  if(pval_approach){
    rownames(arl_track) <- as.character(cl_seq)
    
  }else{
    rownames(arl_track) <- as.character(c(upper.bound.ucl,cl_seq))
    
  }

  return(arl_track)
}

EWMA_Find_CL_AD <- function(ic_data = rnorm(500),
                             lambda = .05,
                             control_limit = .06,
                             m = 5,
                             ICdist = "rnorm",
                             IC_dist_ops = NULL,
                             bootstrap_samples = 3000,
                             track_candidates=NULL) {
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_AD(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  i <- 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  len_com <- len_ic +m
  ind <- 1:(len_com-1)
  lenx_t_leny <- len_ic*m
  f_ic <- ecdf(sorted_ic_data)
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    D_n <- Fast_TS_AD(sorted_ic_data=sorted_ic_data, new_data=data, 
                       lenx=len_ic, f_ic=f_ic,
                       i=ind, len_com=len_com, lenx_t_leny=lenx_t_leny)
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


Fast_TS_CVM <- function(sorted_ic_data, new_data, lenx, i, j, leny){
  # Internal use ONLY 
  ranks <- rank(c(sorted_ic_data, new_data))
  U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
  STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
}

Fast_TS_AD <- function(sorted_ic_data, new_data, lenx, i, len_com, lenx_t_leny, f_ic){
  # form new_data cdf
  M <- lenx * f_ic(sort.int(c(sorted_ic_data, new_data))) # sort before?
  M <- M[-len_com]
  STAT <- 1/(lenx_t_leny)*sum((M * len_com -lenx *i)^2/(i*(len_com-i)))
}

Slow_AD_Two_Sample = function(x,y){
  m = length(y)
  n = length(x)
  N = m + n
  sortx = sort(x)
  combined.sample.sorted = sort(c(x,y))
  A_n = 0
  M = numeric(N-1)
  for(i in 1:(N-1)){
    M = sum(sortx <=combined.sample.sorted[i])
    A_n = (M*N - n*i)^2/(i*(N-i)) + A_n
  }
  M
  A_n = A_n/(m*n)
  #print(A_n)
  A_n

}


Fast_TS_AD_Rev = function(sorted_ic_data, y, m, n, N ){
  # m = length(y)
  # n = length(x)
  # N = m + n
  # sortx = sort(x)
  combined.sample.sorted = sort(c(sorted_ic_data,y))
  A_n = 0
  M = numeric(N-1)
  for(i in 1:(N-1)){
    M = sum(sorted_ic_data <=combined.sample.sorted[i]) # potentially doing way more than needed
    A_n = (M*N - n*i)^2/(i*(N-i)) + A_n
  }
  M
  A_n = A_n/(m*n)
  #print(A_n)
  A_n
  
}

Find_CL_AD_Window <- function(ic_data = rnorm(500),
                                   lambda = 30,
                                   control_limit = .06,
                                   m = 5,
                                   ICdist = "rnorm",
                                   IC_dist_ops = NULL,
                                   method="constant_ad",
                                   track_candidates=NULL) {
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  # Mean of A_mn = 1
  chosen_window=lambda
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)
  ind <- 1:len_ic
  SD_A2= sqrt((2/3)*(pi^2 - 9))
  
  # bootstrap for initial 3 variances
  # var_correct_samp = 3
  # ad.sd = numeric(3)
  # 
  # for(i in 1:var_correct_samp){
  #   ad.boot = numeric(200)
  #   for(j in 1:1000){
  #     b.samp = sample(sorted_ic_data, m*i)
  #     ad.boot[j] =Fast_TS_AD_Rev(sorted_ic_data=sorted_ic_data, y=b.samp, 
  #                                m = m*i, 
  #                                n = len_ic, 
  #                                N = len_ic + m*i)
  #   }
  #   ad.sd[i] = sd(ad.boot)
  # }
  # print(ad.sd)
  i <- 1
  
  while (tail(u, 1) < control_limit) {
    
    new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    data <- cbind(data, new_data) # master data
    num_active_batches <- dim(data)[2]
    amt_data = num_active_batches*m
    
    STAT <- Fast_TS_AD_Rev(sorted_ic_data=sorted_ic_data, y=data, 
                           m = amt_data, 
                           n = len_ic, N =len_ic + amt_data) 
    #print(STAT)
    #print("hello")
    if(amt_data <= 20){
      STAT.Standard = STAT / SD_A2   # ad.sd[i]
    }else{
      STAT.Standard = STAT/ SD_A2
    }
    
    u_nK <- STAT.Standard - 1  # 1 = mean of A_mn
    
    # lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    
    # num_keep <- chosen_window
    # store_keep <- append(store_keep, num_keep)
    # store_remove <- append(store_remove, num_remove)
    # print(paste("On Batch", i))
    # print(paste("Number removed equals ", num_remove))
    # print(paste("Amount of batches used equals ", num_keep))
    if(num_active_batches > chosen_window){
      data = data[,-1]
    }
    print(num_active_batches)
    i <- i + 1
  }
  
  i = i - 1
  
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


Find_CL_AD_Window_OOC <- function(ic_data = rnorm(500),
                                   lambda = 30,
                                   control_limit = .06,
                                   m = 5,
                                   ICdist = "rnorm", IC_dist_ops = NULL,
                                   OOCdist="rnorm", OOC_dist_ops=NULL, tau=0,
                                  method="constant_ad"
                                   ) {
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  # Mean of A_mn = 1
  chosen_window=lambda
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)
  ind <- 1:len_ic
  SD_A2= sqrt((2/3)*(pi^2 - 9))
  
  # bootstrap for initial 3 variances
  # var_correct_samp = 3
  # ad.sd = numeric(3)
  # 
  # for(i in 1:var_correct_samp){
  #   ad.boot = numeric(200)
  #   for(j in 1:1000){
  #     b.samp = sample(sorted_ic_data, m*i)
  #     ad.boot[j] =Fast_TS_AD_Rev(sorted_ic_data=sorted_ic_data, y=b.samp, 
  #                                m = m*i, 
  #                                n = len_ic, 
  #                                N = len_ic + m*i)
  #   }
  #   ad.sd[i] = sd(ad.boot)
  # }
  # print(ad.sd)
  i <- 1
  
  while (tail(u, 1) < control_limit) {
    if(i <= tau){
      new_data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 
      
    }else{
      new_data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    data <- cbind(data, new_data) # master data
    num_active_batches <- dim(data)[2]
    amt_data = num_active_batches*m
    
    STAT <- Fast_TS_AD_Rev(sorted_ic_data=sorted_ic_data, y=data, 
                           m = amt_data, 
                           n = len_ic, N =len_ic + amt_data) 
    #print(STAT)
    #print("hello")
    if(amt_data <= 20){
      STAT.Standard = STAT / SD_A2   # ad.sd[i]
    }else{
      STAT.Standard = STAT/ SD_A2
    }
    
    u_nK <- STAT.Standard - 1  # 1 = mean of A_mn
    
    # lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    
    # num_keep <- chosen_window
    # store_keep <- append(store_keep, num_keep)
    # store_remove <- append(store_remove, num_remove)
    # print(paste("On Batch", i))
    # print(paste("Number removed equals ", num_remove))
    # print(paste("Amount of batches used equals ", num_keep))
    if(num_active_batches > chosen_window){
      data = data[,-1]
    }
    print(num_active_batches)
    i <- i + 1
  }
  
  
  i = i - 1
  # # If we have a large upper bound, we can calculate all our upper bounds below that number
  # if(!is.null(track_candidates)){
  #   track_ucl <-numeric(length(track_candidates))
  #   
  #   for(j in seq_along(track_candidates)){
  #     track_ucl[j] <- as.numeric(which.min(track_candidates[j] > u)) - 1 #investigate later
  #   }
  #   names(track_ucl) <- as.character(track_candidates)
  #   return(list("Time OOC"=i, "Lower CLs"=track_ucl))
  #   
  # }
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

EWMA_Find_CL_KS <- function(ic_data = rnorm(500),
                             lambda = .05,
                             control_limit = .06,
                             m = 5,
                             ICdist = "rnorm",
                             IC_dist_ops = NULL,
                             bootstrap_samples = 3000,
                             track_candidates=NULL) {
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_KS(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  i <- 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  ## sorted_ic_data <- sort(ic_data) # No longer needed for KS
  ic_fun <- ecdf(ic_data)
  len_ic <- length(ic_data)
  
  
  ind <- 1:len_ic
  
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    D_n <- Fast_TS_KS(ic_fun=ic_fun, new_data=data,
                       j=j, m=m)
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
  lenx_t_leny <- lenx*leny
  lenx_p_leny <- lenx+leny
  
  for (k in 1:bootstrap_samples) {
    sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
    ranks <- rank(c(sort_ic, sampled_data))
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):lenx_p_leny] - j)^2)
    D_n[k] <- U / ((lenx_t_leny)*(lenx_p_leny)) - (4 * lenx_t_leny - 1)/(6 * (lenx_p_leny))
  }
  d0 = mean(D_n)
}




Fast_Bootstrap_KS <- function(ic_data, m, bootstrap_samples){
  # Tested 9/15/2016
  fhat_ic <- ecdf(ic_data)
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    D_n[i] <-  max(f0 - (j - 1) / m, j / m - f0)
  }
  d0=mean(D_n)
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
  for(k in 1:bootstrap_samples) {
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
  # This is termed "experimental due to actually not simulating from the in-control distribution
  D_n <- numeric(bootstrap_samples)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  i <- 1:lenx
  leny <- m
  j <- 1:leny

  # Newly simulated data
  rIC <- get(ICdist, mode = "function", envir = parent.frame())


  for(k in 1:bootstrap_samples) {
    sampled_data <- do.call(rIC, c(list(m), ic_params))  #NOW A MISNOMER
    ranks <- rank(c(sort_ic, sampled_data))
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    D_n[k] <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
  }
  D_n
}

Fast_Bootstrap_AD <- function(ic_data, m, bootstrap_samples){
  # This is inaccurate (5/5/2016)
  # Possibly because of ties
  D_n <- numeric(bootstrap_samples)
  f_ic <- ecdf(ic_data)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  # i <- 1:lenx
  leny <- m
  # j <- 1:leny
  len_com <- lenx + leny
  i <- 1:(len_com - 1)
  lenx_t_leny <- lenx * leny
  for (k in 1:bootstrap_samples) {
    sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
    M <- lenx * f_ic(sort.int(c(sort_ic, sampled_data))) # sort before?
    M <- M[-len_com]
    D_n[k] <- 1/(lenx_t_leny)*sum((M * len_com -lenx *i)^2/(i*(len_com-i)))
  }
  d0 = mean(D_n)
}

In_Control_Dist_AD_Exp <- function(ic_data, m, bootstrap_samples, ICdist, ic_params){
  # This is termed "experimental due to actually not simulating from the in-control distribution
  D_n <- numeric(bootstrap_samples)
  f_ic <- ecdf(ic_data)
  sort_ic <- sort(ic_data)
  lenx <- length(ic_data)
  # i <- 1:lenx
  leny <- m
  # j <- 1:leny
  len_com <- lenx + leny
  i <- 1:(len_com - 1)
  lenx_t_leny <- lenx*leny
  
  # Newly simulated data
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  for(k in 1:bootstrap_samples) {
    sampled_data <- do.call(rIC, c(list(m), ic_params))  #NOW A MISNOMER
    M <- lenx * f_ic(sort.int(c(sort_ic, sampled_data))) # sort before?
    M <- M[-len_com]
    D_n[k] <- 1/(lenx_t_leny)*sum((M * len_com -lenx *i)^2/(i*(len_com-i)))
  }
  D_n
}


In_Control_Dist_KS_Exp <- function(ic_data, m, bootstrap_samples, ICdist, ic_params){
  # 
  D_n <- numeric(bootstrap_samples)
  fhat_ic <- ecdf(ic_data)
  j <- 1:m
  # Newly simulated data
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  for(k in 1:bootstrap_samples) {
    sampled_data <- do.call(rIC, c(list(m), ic_params))  #NOW A MISNOMER
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    D_n[k] <-  max(f0 - (j - 1) / m, j / m - f0)
  }
  D_n
}




Fast_TS_KS <- function(ic_fun, new_data, j, m){
 # Internal use ONLY 
  new_data <- sort(new_data)
  f0 <- ic_fun(new_data)
  max(f0-(j-1)/m, j/m- f0)
}




EWMA_CVM_OOC <- function(ic_data=rnorm(500), lambda=.05, control_limit, m=5, exact=FALSE, 
                          tau=0, 
                          ICdist="rnorm", IC_dist_ops=NULL,
                          OOCdist="rnorm", OOC_dist_ops=NULL,
                          bootstrap_samples=3000){
  # TODO: Add multicore support
  
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



EWMA_AD_OOC <- function(ic_data=rnorm(500), lambda=.05, control_limit, m=5, exact=FALSE, 
                         tau=0, 
                         ICdist="rnorm", IC_dist_ops=NULL,
                         OOCdist="rnorm", OOC_dist_ops=NULL,
                         bootstrap_samples=3000){
  # TODO: Add multicore support
  
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_AD(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
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
  len_com <- len_ic + m
  ind <- 1:(len_com-1)
  lenx_t_leny <- len_ic*m
  f_ic <- ecdf(sorted_ic_data)
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 
      
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    # Calculate new value of test statisic
    D_n <- Fast_TS_AD(sorted_ic_data=sorted_ic_data, new_data=data, 
                      lenx=len_ic, 
                      i=ind, len_com=len_com, lenx_t_leny=lenx_t_leny, f_ic=f_ic)
    # print(D_n)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
  }
  return(list("control_stats" = u, "Time OOC" = i))
}

EWMA_KS_OOC <- function(ic_data=rnorm(500), lambda=.05, control_limit, m=5, exact=FALSE, 
                         tau=0, 
                         ICdist="rnorm", IC_dist_ops=NULL,
                         OOCdist="rnorm", OOC_dist_ops=NULL,
                         bootstrap_samples=3000){
  # TODO: Add multicore support
  
  # Calculate d0 via bootstrapping
  d0 <- Fast_Bootstrap_KS(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  ic_fun <- ecdf(ic_data)
  len_ic <- length(ic_data)

  
  while (tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 
      
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    # Calculate new value of test statisic
    D_n <- Fast_TS_KS(ic_fun=ic_fun, new_data=data,
                      j=j, m=m)
    # print(D_n)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
  }
  return(list("control_stats" = u, "Time OOC" = i))
}

##### Core Functions End Here #####




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
