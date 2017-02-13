# Pruning Functions

find_cvm_pval_asym <- function(n, STAT){
  mu <- 1/6
  lambda <- 5*n /(24*n-18)
  STAT<- 1-pinvgauss(STAT, mu, lambda)
}

find_ad_pval_asym <- function(STAT){
  # from Marsalia paper
  if(STAT <=0) return("ERROR, STAT out of bounds")
  if(STAT<2 & STAT > 0){
   pval <-  1- STAT^(-.5)*exp(-1.2337141/STAT)*(2.00012 + (0.247105 -(.0649821 -(.0347962 -(.0116720 -.00168691*STAT)*STAT)*STAT)*STAT)*STAT)
  }else{
   pval <- 1- exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056-.0003146*STAT)*STAT)*STAT)*STAT)*STAT)) 
  }
  # mu <- 1/6
  # lambda <- 5*n /(24*n-18)
  # STAT<- 1-pinvgauss(STAT, mu, lambda)
}


EWMA_Find_CL_CVM_Pval <- function(ic_data = rnorm(500),
                                  lambda = .05,
                                  control_limit = .01,
                                  m = 5,
                                  ICdist = "rnorm",
                                  IC_dist_ops = NULL,
                                  track_candidates=NULL,
                                  diag=FALSE,
                                  method="partha",
                                  run_until=FALSE){
  # run_until: diagnostic option. It should always be set to false
  # diag: diagnostic option
  # track_candidates: diagnostic option. Should be true/false, but it isn't. oh well
  
  require(statmod, quietly=TRUE) # needed to find cvm pvalue approximation

  constant <- lambda # Makes things more compatible later
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  f_ic <- ecdf(sorted_ic_data)

  store_keep <- NULL
  store_remove <- NULL
  
  # Pruning methods put here. If outside function won't compile right
  source("functions_prune_internal.R", local=TRUE)
  if(run_until){
    control_limit <- 0
  }
  
  while (tail(u, 1) > control_limit) {
    # Calculate new value of test statisic
    new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    data <- cbind(data, new_data) # master data
    num_active_batches <- dim(data)[2]
    #if(dim(data)[2] > 2000) stop("data has grown very large")
    # Pruning Function Here
    if(i==1){
      STAT <- cvm.res(x=data, y=f_ic)
      # find p-value 
      pval <- find_cvm_pval_asym(n=m, STAT=STAT)
      if(method=="ewma"){
        ewma_stat <- pval  # average p-val?
      }
    }else{
      # on the previous steps data + new data, run a cvm test
      if(method=="lookback"){
        STAT.big <- cvm.res(x=data, y=f_ic)
        pval.big <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT.big)
        STAT.small <- cvm.res(x=data[,-1], y=f_ic)
        pval.small <- find_cvm_pval_asym(n=m*(num_active_batches-1), STAT=STAT.small)
      }else {
        STAT <- cvm.res(x=data, y=f_ic)
        pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      }
      # find amount of previous data to remove
      
      if(method=="lookback"){
        num_remove <- prune(pval=pval.big, method=method, pval.small=pval.small)
        pval = pval.big #idk, maybe
      }else{
        num_remove <- prune(pval=pval, method=method)
      }
      
      # For an ewma p-value method, what we gon do?
      if(method=="ewma"){
        ewma_stat = constant * pval + (1 - constant) * u_nK
      }
      
      # for diagnostics
      num_keep <- num_active_batches - num_remove
      store_keep <- append(store_keep, num_keep)
      store_remove <- append(store_remove, num_remove)
      # print(paste("On Batch", i))
      # print(paste("Number removed equals ", num_remove))
      # print(paste("Amount of batches used equals ", num_keep))
      if(num_remove > 0){
        data <- data[,-1:-num_remove]
      }
    }
    # update control statistic 
    if(method=="ewma"){
      u_nK = ewma_stat
    }else{    
      u_nK <- pval
    }
    u <- append(u, u_nK)
    i <- i + 1

    if(run_until){
      if(i == 300) u[length(u)] <- -.2
    }
  }
  
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  
  # For the p-value approach, we can calclulate for the SMALLEST p-value
  if(!is.null(track_candidates)){
    track_ucl <-numeric(length(track_candidates))
    
    for(j in seq_along(track_candidates)){
      
      #print(as.numeric(which.min(track_candidates[j] < u)))
      track_ucl[j] <- as.numeric(which.min(track_candidates[j] < u)) - 1 #investigate later
    }
    names(track_ucl) <- as.character(track_candidates)
    return(list("Time OOC"=i, "Lower CLs"=track_ucl))
  }
  if(diag){
    return(list("control_stats" = u, "Time OOC" = i, "Track Removed"=store_remove,
                "Track Keep"=store_keep))
  }
  return(list("control_stats" = u, "Time OOC" = i))
}

EWMA_Find_CL_CVM_Pval_OOC <- function(ic_data=rnorm(500), 
                                      lambda=.05, 
                                      control_limit=.01, 
                                      m=5, 
                                      tau=0, 
                                      ICdist="rnorm", IC_dist_ops=NULL,
                                      OOCdist="rnorm", OOC_dist_ops=NULL,
                                      method="partha"){

  # 11/1/2016 -- Checked for accuracy
  require(statmod, quietly=TRUE) # Needed for asymptotic p-value
  constant <- lambda # Makes things more compatible later
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  f_ic <- ecdf(sorted_ic_data)
  
  # For diagnostic functions
  store_keep <- NULL
  store_remove <- NULL
  
  # pruning function
  # Pruning methods put here. If outside function won't compile right
  source("functions_prune_internal.R", local=TRUE)
  

  while (tail(u, 1) > control_limit) {
    if(i <= tau){
      new_data <- do.call(rIC, c(n=list(m), IC_dist_ops))
    }else{
      # example calls:
      # rOOC= standard_chi, list(df=1, center_mean=1, center_sd=1)
      new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
    }
    data <- cbind(data, new_data) # master data
    num_active_batches <- dim(data)[2]
    if(i==1){
      STAT <- cvm.res(x=data, y=f_ic)
      # find p-value 
      pval <- find_cvm_pval_asym(n=m, STAT=STAT)
      if(method=="ewma"){
        ewma_stat <- pval  # Is this right? I feel like this is going to bias results!
      }
    }else{
      # on the previous steps data + new data, run a cvm test
      if(method=="lookback"){
        STAT.big <- cvm.res(x=data, y=f_ic)
        pval.big <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT.big)
        STAT.small <- cvm.res(x=data[,-1], y=f_ic)
        pval.small <- find_cvm_pval_asym(n=m*(num_active_batches-1), STAT=STAT.small)
      }else {
        STAT <- cvm.res(x=data, y=f_ic)
        pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      }
      # find amount of previous data to remove
      
      if(method=="lookback"){
        num_remove <- prune(pval=pval.big, method=method, pval.small=pval.small)
        pval = pval.big #idk, maybe
      }else{
        num_remove <- prune(pval=pval, method=method)
      }
      
      # For an ewma p-value method, what we gon do?
      if(method=="ewma"){
        ewma_stat = constant * pval + (1 - constant) * u_nK
      }
      
      # for diagnostics
      num_keep <- num_active_batches - num_remove
      store_keep <- append(store_keep, num_keep)
      store_remove <- append(store_remove, num_remove)
      # print(paste("On Batch", i))
      # print(paste("Number removed equals ", num_remove))
      # print(paste("Amount of batches used equals ", num_keep))
      if(num_remove > 0){
        data <- data[,-1:-num_remove]
      }
    }
    # update control statistic 
    if(method=="ewma"){
      u_nK = ewma_stat
    }else{    
      u_nK <- pval
    }
    u <- append(u, u_nK)
    i <- i + 1

  }
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
}


EWMA_Find_CL_CVM_Prune <- function(ic_data = rnorm(500),
                                   lambda = .05,
                                   control_limit = .06,
                                   m = 5,
                                   ICdist = "rnorm",
                                   IC_dist_ops = NULL,
                                   bootstrap_samples = 3000,
                                   track_candidates=NULL,
                                   constant=.5,
                                   pval_table){
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  num_history <- 15
  tr <- 1:num_history
  weights <- rep(.9^tr, each=num_history)
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  how_far_back <- 1
  
  while (tail(u, 1) < control_limit) {
    # if(i <= num_history){
    #   tr_temp <- 1:i
    #   new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    #   data <- cbind(data, new_data)
    #   if(i==1){
    #     data <- data[, -1]
    #   }
    #   weights <- rep(.9^tr_temp, each=m)
    #   len.com.temp <- len_ic + dim(data)[1]*dim(data)[2]
    #   ind.temp <- 1:(len_com-1)
    #   lenx.t.leny.temp <- len_ic * dim(data)[1]*dim(data)[2]
    # }
    # if(i <= num_history){
    #   new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    #   data <- data[, -1]
    #   data <- cbind(data, new_data)
    # }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    data <- cbind(data, new_data) # master data
    if(dim(data)[2] > 2000) stop("data has grown very large")
    paste("at step ", i, " the dimension of data is", dim(data)[2])
    # Pruning Function Here
    if(i==1){
      p1 <- 1
      STAT <- cvm.res(x=data, y=f_ic)
      print("i=1")
      # lookup p-value from table
      pval <- find_pval(N=as.character(m), STAT=STAT)
      # Error catch if N is not found in table needs to go HERE
    }else{
      
      # on the previous steps data + new data, run a cvm test
      STAT <- cvm.res(x=data, y=f_ic)
      pval <- find_pval(N=as.character(m), STAT=STAT, pval_table=pval_table)
      # re-adjust effective (active) data
      # obviously, clean this up later to make it look nicer
      percent_prune <- min(.2, ((pval-constant)/(1-constant))^2)
      print(((pval-constant)/(1-constant))^2)
      percent_prune
      percent_keep <- 1 - percent_prune
      num_keep <- ceiling(dim(data)[2] * percent_keep)
      num_remove <- dim(data[2])-num_keep
      if(num_remove>0){
        data <- data[,-1:num_remove]
      }
    }
    
    #D_n <- cvm.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- lambda * (STAT - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    if(i > 20){
      u[length(u)] <- control_limit + 1
    }
  }
  i <- i -1 # this adjust for the time out of control!
  
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


EWMA_Find_CL_AD_Weight <- function(ic_data = rnorm(500),
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
  u <- 0 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  num_history <- 15
  tr <- 1:num_history
  weights <- rep(.9^tr, each=num_history)
  
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  len_com <- len_ic + m*num_history
  ind <- 1:(len_com-1)
  lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  while (tail(u, 1) < control_limit) {
    if(i <= num_history){
      tr_temp <- 1:i
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- cbind(data, new_data)
      if(i==1){
        data <- data[, -1]
      }
      weights <- rep(.9^tr_temp, each=m)
      len.com.temp <- len_ic + dim(data)[1]*dim(data)[2]
      ind.temp <- 1:(len_com-1)
      lenx.t.leny.temp <- len_ic * dim(data)[1]*dim(data)[2]
    }
    if(i <= num_history){
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- data[, -1]
      data <- cbind(data, new_data)
    }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    D_n <- ad.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- lambda * (D_n - d0) + (1 - lambda) * tail(u, 1)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    if(i > 20){
      u[length(u)] <- control_limit + 1
    }
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

# Pruning Functions

