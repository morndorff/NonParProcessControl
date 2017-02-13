# Weighted CVM w/ weighted ecdf

Find_CL_CVM_Weight2 <- function(ic_data = rnorm(500),
                                    lambda = .05,
                                    control_limit = .06,
                                    m = 5,
                                    ICdist = "rnorm",
                                    IC_dist_ops = NULL,
                                    track_candidates=NULL, 
                                    W=30,
                                    method="weighted") {
  # See below for track_candidates explanation
  require(spatstat)
  require(statmod)
  
  if(method!="weighted") stop("Improper method specified!")
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  #num_history <- 15
  #tr <- 1:num_history
  #weights <- rep(.9^tr, each=num_history)
  
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  while(tail(u, 1) > control_limit) {
    if(i <= W){
      # calculate weights here
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      #print(data)
      #print(i)
      #print(current.weight)
      #print(data)
      STAT <- cvm.res.weight.internal(x=data, y=f_ic, w=current.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }else if(i >= W){
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      # print(i)
      # print(W)
      # print(current.weight)
      # print(dim(data)[2])
      STAT <- cvm.res.weight.internal(x=data, y=f_ic, w=maxed.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    # D_n <- cvm.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- D_n
    print(D_n)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    # if(i > 40){
    #   u[length(u)] <- control_limit - 1
    # }
  }
  
  i = i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  
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
  return(list("control_stats" = u, "Time OOC" = i))
}



Find_CL_CVM_Weight2_OOC <- function(ic_data=rnorm(500), 
                                    lambda=.05, 
                                    control_limit=.01, 
                                    m=5, 
                                    tau=0, 
                                    ICdist="rnorm", IC_dist_ops=NULL,
                                    OOCdist="rnorm", OOC_dist_ops=NULL,
                                    W=30,
                                    method="weighted") {
  # See below for track_candidates explanation
  require(spatstat)
  require(statmod)
  
  if(method!="weighted") stop("Improper method specified!")
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  #num_history <- 15
  #tr <- 1:num_history
  #weights <- rep(.9^tr, each=num_history)
  
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  while(tail(u, 1) > control_limit) {
    if(i <= W){
      # calculate weights here
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      if(i <=tau){
        new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
      }
      

      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      #print(data)
      #print(i)
      #print(current.weight)
      #print(data)
      STAT <- cvm.res.weight.internal(x=data, y=f_ic, w=current.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }else if(i >= W){
      if(i <=tau){
        new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
      }
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      # print(i)
      # print(W)
      # print(current.weight)
      # print(dim(data)[2])
      STAT <- cvm.res.weight.internal(x=data, y=f_ic, w=maxed.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    # D_n <- cvm.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- D_n
    # print(D_n)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    # if(i > 40){
    #   u[length(u)] <- control_limit - 1
    # }
  }
  
  i = i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
  
}




Find_CL_CVM_Weight_Kernel <- function(ic_data = rnorm(500),
                                lambda = .05,
                                control_limit = .06,
                                m = 5,
                                ICdist = "rnorm",
                                IC_dist_ops = NULL,
                                track_candidates=NULL, 
                                W=30,
                                method="weighted_kernel") {
  # See below for track_candidates explanation
  require(spatstat)
  require(statmod)
  
  if(method!="weighted_kernel") stop("Improper method specified!")
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  #num_history <- 15
  #tr <- 1:num_history
  #weights <- rep(.9^tr, each=num_history)
  
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  while(tail(u, 1) > control_limit) {
    if(i <= W){
      # calculate weights here
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      #print(data)
      #print(i)
      #print(current.weight)
      #print(data)
      STAT <- cvm.res.weight.internal.kernel(x=data, y=f_ic, w=current.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }else if(i >= W){
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      # print(i)
      # print(W)
      # print(current.weight)
      # print(dim(data)[2])
      STAT <- cvm.res.weight.internal.kernel(x=data, y=f_ic, w=maxed.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    # D_n <- cvm.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- D_n
    print(D_n)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    # if(i > 40){
    #   u[length(u)] <- control_limit - 1
    # }
  }
  
  i = i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  
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
  return(list("control_stats" = u, "Time OOC" = i))
}



Find_CL_CVM_Weight_Kernel_OOC <- function(ic_data=rnorm(500), 
                                    lambda=.05, 
                                    control_limit=.01, 
                                    m=5, 
                                    tau=0, 
                                    ICdist="rnorm", IC_dist_ops=NULL,
                                    OOCdist="rnorm", OOC_dist_ops=NULL,
                                    W=30,
                                    method="weighted_kernel") {
  # See below for track_candidates explanation
  require(spatstat)
  require(statmod)
  
  if(method!="weighted_kernel") stop("Improper method specified!")
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 1 # initialize first value at 1 (p-value approach)
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  #num_history <- 15
  #tr <- 1:num_history
  #weights <- rep(.9^tr, each=num_history)
  
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  while(tail(u, 1) > control_limit) {
    if(i <= W){
      # calculate weights here
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      if(i <=tau){
        new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
      }
      
      
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      #print(data)
      #print(i)
      #print(current.weight)
      #print(data)
      STAT <- cvm.res.weight.internal.kernel(x=data, y=f_ic, w=current.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }else if(i >= W){
      if(i <=tau){
        new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
      }
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      # print(i)
      # print(W)
      # print(current.weight)
      # print(dim(data)[2])
      STAT <- cvm.res.weight.internal.kernel(x=data, y=f_ic, w=maxed.weight)
      D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    # data <- do.call(rIC, c(list(m), IC_dist_ops))
    # Calculate new value of test statisic
    # D_n <- cvm.weighted(x=data, y=f_ic, w=weights)
    # print(D_n)
    u_nK <- D_n
    # print(D_n)
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    # if(i > 40){
    #   u[length(u)] <- control_limit - 1
    # }
  }
  
  i = i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
  
}
# Weighted CVM w/ weighted ecdf

EWMA_Find_CL_CVM_Weight <- function(ic_data = rnorm(500),
                                    lambda = .05,
                                    W = 30,
                                    control_limit = .06,
                                    m = 5,
                                    ICdist = "rnorm",
                                    IC_dist_ops = NULL,
                                    bootstrap_samples = 10000,
                                    track_candidates=NULL,
                                    method="weighted_bs",
                                    cvm_approx= f_mem,
                                    diag=FALSE) {
  #cvm_approx = memoised function! 
  require(spatstat)
  require(memoise)
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
  # num_history <- 15
  # tr <- 1:num_history
  # weights <- rep(.9^tr, each=num_history)
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  
  while (tail(u, 1) < control_limit) {
    new_data <- do.call(rIC, c(list(m), IC_dist_ops))
    
    if(i<=W){
      data <- cbind(data, new_data) # master data
      num_active_batches <- dim(data)[2]
      
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      # new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      # data <- cbind(data, new_data)
      # num_active_batches = dim(data)[2]
      #if(is.memoised(cvm_approx)) print("we memed")
      #else print("we didn't meme")
      # calculate bootstraped value
      #print(paste("ic_sample_size", len_ic) )
      #print(paste("m", m * num_active_batches))
      #print(paste("bootstrap_samples", bootstrap_samples))
      # print(paste("w" , current.weight))
      # print(environment(cvm_approx))
      ti=proc.time()
      mean_var_approx = cvm_approx(ic_sample_size = len_ic, 
                      m = m * num_active_batches, 
                      bootstrap_samples=bootstrap_samples, 
                      w=current.weight)
      ti2 = proc.time() - ti
      print(ti2)
      d_mean=mean_var_approx[["mean"]]
      d_sd = mean_var_approx[["sd"]]
      D_n <- cvm.res.weight.internal(x=data, y=f_ic, w=current.weight)
      # D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
      
      # cvm.res.weight.internal(x=data, y=f_ic, w=w)
      # STAT <- cvm.res(x=data, y=f_ic)
      # find p-value 
      # pval <- find_cvm_pval_asym(n=m, STAT=STAT)
    }else{
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]

      mean_var_approx = cvm_approx(ic_sample_size = len_ic, 
                      m = m * num_active_batches, 
                      bootstrap_samples=bootstrap_samples, 
                      w=maxed.weight)
      d_mean = mean_var_approx[["mean"]]
      d_sd = mean_var_approx[["sd"]]
      
      D_n <- cvm.res.weight.internal(x=data, y=f_ic, w=maxed.weight)
      # D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    
    u_nK <-  (D_n - d_mean)/d_sd # standardized control statistic
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    if(diag){
    # stop after 20 regardless
      if(i > 50){
       u[length(u)] <- control_limit + 1
      }
    }
  }
  
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  
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



EWMA_Find_CL_CVM_Weight_OOC <- function(ic_data = rnorm(500),
                                    lambda = .05,
                                    W = 30,
                                    control_limit = .06,
                                    m = 5,
                                    ICdist = "rnorm", IC_dist_ops = NULL,
                                    OOCdist="rnorm", OOC_dist_ops=NULL,
                                    bootstrap_samples = 10000,
                                    track_candidates=NULL,
                                    method="weighted_bs", tau=0,
                                    cvm_approx) {
  #cvm_approx = memoised function! 
  require(spatstat)
  require(memoise)
  # See below for track_candidates explanation
  
  # Calculate d0 via bootstrapping
  # d0 <- Fast_Bootstrap_CVM(ic_data=ic_data, m=m, bootstrap_samples=bootstrap_samples)
  
  # Initializing Variables
  data <- NULL
  u <- 0 # initialize first value at 0
  i <- 1
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  
  # Will be parameterized later
  # num_history <- 15
  # tr <- 1:num_history
  # weights <- rep(.9^tr, each=num_history)
  
  # Pre Calculations for functions
  sorted_ic_data <- sort(ic_data)
  len_ic <- length(sorted_ic_data)  
  # len_com <- len_ic + m*num_history
  # ind <- 1:(len_com-1)
  # lenx_t_leny <- len_ic *  m*num_history
  f_ic <- ecdf(sorted_ic_data)
  
  # Precalculations for weights
  maxed.weight = rep((1-lambda)^(W - 1:W), each=m)
  maxed.weight = maxed.weight/sum(maxed.weight)
  
  
  
  while (tail(u, 1) < control_limit) {
    if(i <= tau){
      new_data <- do.call(rIC, c(n=list(m), IC_dist_ops))
    }else{
      # example calls:
      # rOOC= standard_chi, list(df=1, center_mean=1, center_sd=1)
      new_data <- do.call(rOOC, c(n=list(m), OOC_dist_ops))
    }
    if(i<=W){
      data <- cbind(data, new_data) # master data
      num_active_batches <- dim(data)[2]
      
      current.weight = rep((1-lambda)^(i - 1:i), each=m)
      current.weight = current.weight/sum(current.weight)
      # tr_temp <- 1:i
      # new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      # data <- cbind(data, new_data)
      # num_active_batches = dim(data)[2]
      # calculate bootstraped value
      mean_var_approx = cvm_approx(ic_sample_size = len_ic, 
                                   m = m * num_active_batches, 
                                   bootstrap_samples=bootstrap_samples, 
                                   w=current.weight)
      d_mean=mean_var_approx[["mean"]]
      d_sd = mean_var_approx[["sd"]]
      
      D_n <- cvm.res.weight.internal(x=data, y=f_ic, w=current.weight)
      # D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
      
      # cvm.res.weight.internal(x=data, y=f_ic, w=w)
      # STAT <- cvm.res(x=data, y=f_ic)
      # find p-value 
      # pval <- find_cvm_pval_asym(n=m, STAT=STAT)
    }else{
      new_data <- do.call(rIC, c(list(m), IC_dist_ops))
      data <- data[, -1]
      data <- cbind(data, new_data)
      num_active_batches = dim(data)[2]
      
      mean_var_approx = cvm_approx(ic_sample_size = len_ic, 
                                   m = m * num_active_batches, 
                                   bootstrap_samples=bootstrap_samples, 
                                   w=maxed.weight)
      d_mean=mean_var_approx[["mean"]]
      d_sd = mean_var_approx[["sd"]]
      D_n <- cvm.res.weight.internal(x=data, y=f_ic, w=maxed.weight)
      # D_n = find_cvm_pval_asym(n = m* num_active_batches, STAT=STAT) 
    }
    
    u_nK <-  (D_n - d_mean)/d_sd # standardized control statistic
    # print(paste("u=",u_nK))
    u <- append(u, u_nK)
    i <- i + 1
    # stop after 20 regardless
    # if(i > 20){
    #   u[length(u)] <- control_limit + 1
    # }
  }
  
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
}










# The bootstrap, weighted

Fast_Bootstrap_CVM_Weighted <- function(ic_sample_size, m, bootstrap_samples=10000, w=NULL){
  # This is inaccurate (5/5/2016)
  # Possibly because of ties
  if(is.null(w)){
    w = rep(1/m, m)
    print("equal weighting used")
  }
  D_n <- numeric(bootstrap_samples)
  ic_data = 1:ic_sample_size
  #sort_ic <- sort(ic_data)
  # lenx <- length(ic_data)
  # i <- 1:lenx
  y <- ecdf(ic_data)
  
  # leny <- m
  # j <- 1:leny
  # lenx_t_leny <- lenx*leny
  # lenx_p_leny <- lenx+leny
  
  for (k in 1:bootstrap_samples) {
    sampled_data <- sort.int(sample(ic_data, replace = T, size = m))
    D_n[k] = cvm.res.weight.internal(x= sampled_data, y=y, w=w)
    # ranks <- rank(c(sort_ic, sampled_data))
    # U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):lenx_p_leny] - j)^2)
    # D_n[k] <- U / ((lenx_t_leny)*(lenx_p_leny)) - (4 * lenx_t_leny - 1)/(6 * (lenx_p_leny))
  }
  d0 = mean(D_n)
  d_sd = sd(D_n)
  return(list("mean"=d0, "sd"=d_sd))
}



cvm.res.weight.internal.kernel <- function(x, y, w, ...){
  # want to recieve weights from outside the function, since it will be easier
  
  
  # x: numeric vector
  # y: numeric vector or prob distribution e.g. "pnorm" 
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    x <- sort(x)
    y <- sort(y)
    leny <- length(y)
    # com <- c(x,y)
    # ranks <- rank(com)
    # lenx <- length(x)
    # leny <- length(y)
    # i <- 1:lenx
    # j <- 1:leny
    # # diagnostic 
    # print((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # plot((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
    # print(((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx)))
    f_ic <- ecdf(x)
    f_new <- ecdf(y)
    change <- sum((f_ic(y)-f_new(y))^2) + sum((f_ic(x)-f_new(x))^2) #lol how do I normalize
    change <- change * (lenx * leny)/(lenx + leny)^2
    return(STAT)
  }
  # One Sample Test  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  F_x <- y(x, ...)
  i <- 1:lenx
  # November update. Need to get weighting right here
  # F_w = ewcdf(x, w) # W MUST BE NORMALIZED
  f_w = density(x, weights=w)
  # get y
  testy = cumsum(f_w$y)/sum(f_w$y)
  F_w = approxfun(f_w$x, testy, yleft=0, yright=1)

  # W also needs taht batch size adjustment to be done beforehand
  # print("weighted cdf")
  # print(F_w(x))
  
  
  STAT <- (1/(12*lenx)) + sum((F_x - (F_w(x) - (.5/lenx) ) )^2)
  # STAT <- (1/(12*lenx)) + sum( (F_x -(1:lenx - .5)/ lenx)^2) old
  STAT
}

cvm.res.weight.internal <- function(x, y, w, ...){
  # want to recieve weights from outside the function, since it will be easier
  
  # if(sum(w)!=1) stop("weights do not sum to 1")
  
  # x: numeric vector
  # y: numeric vector or prob distribution e.g. "pnorm" 
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    x <- sort(x)
    y <- sort(y)
    leny <- length(y)
    # com <- c(x,y)
    # ranks <- rank(com)
    # lenx <- length(x)
    # leny <- length(y)
    # i <- 1:lenx
    # j <- 1:leny
    # # diagnostic 
    # print((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # plot((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    # STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
    # print(((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx)))
    f_ic <- ecdf(x)
    f_new <- ecdf(y)
    change <- sum((f_ic(y)-f_new(y))^2) + sum((f_ic(x)-f_new(x))^2) #lol how do I normalize
    change <- change * (lenx * leny)/(lenx + leny)^2
    return(STAT)
  }
  # One Sample Test  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  F_x <- y(x, ...)
  i <- 1:lenx
  # November update. Need to get weighting right here
  F_w = ewcdf(x, w) # W MUST BE NORMALIZED
  # W also needs taht batch size adjustment to be done beforehand
  # print("weighted cdf")
  # print(F_w(x))
  STAT <- (1/(12*lenx)) + sum((F_x - (F_w(x) - (.5/lenx) ) )^2)
  # STAT <- (1/(12*lenx)) + sum( (F_x -(1:lenx - .5)/ lenx)^2) old
  STAT
}
