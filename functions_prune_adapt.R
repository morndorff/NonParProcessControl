
Find_CL_CVM_Pval_Adapt <- function(ic_data = rnorm(500),
                                   lambda = 30,
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
  
  require(statmod) # needed to find cvm pvalue approximation
  
  W <- lambda # Makes things more compatible later
  
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
    }else{
      # on the previous steps data + new data, run a cvm test
      
      STAT <- cvm.res(x=data, y=f_ic)
      pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      
      # find amount of previous data to remove
      
      # compare this pvalue to the previous iterations pval
      pval.prev = u_nK 
      
      if(i < W){
        if(pval.prev > pval){
          # W = W # this should result in the window growing to size W
          # small=TRUE      
          num_remove <- 1
        }else{
          # W = W + 1 
          # small=FALSE
          num_remove = 0
        }
      }else{
        # In theory, this sets the window to size 30
        if(pval.prev > pval){
          # W = W -1
          # small=TRUE      
          num_remove <- 2
          if(num_active_batches < 2){
           num_remove = 0 # minimum window size
          }
        }else{
          # W = W + 1 
          # small=FALSE
          num_remove = 0
        }
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
    u_nK <- pval
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

Find_CL_CVM_Pval_Adapt_OOC <- function(ic_data=rnorm(500), 
                                       lambda=.05, 
                                       control_limit=.01, 
                                       m=5, 
                                       tau=0, 
                                       ICdist="rnorm", IC_dist_ops=NULL,
                                       OOCdist="rnorm", OOC_dist_ops=NULL,
                                       method="partha"){
  # 11/1/2016 -- Checked for accuracy
  require(statmod) # Needed for asymptotic p-value
  W <- lambda # Makes things more compatible later
  
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
      pval <- find_cvm_pval_asym(n=m, STAT=STAT)
    }else{
      
      # on the previous steps data + new data, run a cvm test
      STAT <- cvm.res(x=data, y=f_ic)
      pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      
      # find amount of previous data to remove
      
      # compare this pvalue to the previous iterations pval
      pval.prev = u_nK 
      
      if(i < W){
        if(pval.prev > pval){
          # W = W # this should result in the window growing to size W
          # small=TRUE      
          num_remove <- 1
        }else{
          # W = W + 1 
          # small=FALSE
          num_remove = 0
        }
      }else{
        # In theory, this sets the window to size 30
        if(pval.prev > pval){
          # W = W -1
          # small=TRUE      
          num_remove <- 2
          if(num_active_batches < 2){
            num_remove = 0 # minimum window size
          }
        }else{
          # W = W + 1 
          # small=FALSE
          num_remove = 0
        }
      }
      
      # for diagnostics
      num_keep <- num_active_batches - num_remove
      store_keep <- append(store_keep, num_keep)
      store_remove <- append(store_remove, num_remove)
      print(paste("On Batch", i))
      print(paste("Number removed equals ", num_remove))
      print(paste("Amount of batches used equals ", num_keep))
      if(num_remove>0){
        data <- data[,-1:-num_remove]
      }
    }
    u_nK <- pval
    u <- append(u, u_nK)
    i <- i + 1
    
  }
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
}


Find_CL_CVM_Pval_DWindow <- function(ic_data = rnorm(500),
                                  lambda = .05,
                                  control_limit = .01,
                                  m = 5,
                                  ICdist = "rnorm",
                                  IC_dist_ops = NULL,
                                  track_candidates=NULL,
                                  diag=FALSE,
                                  method="partha",
                                  run_until=FALSE,
                                  W_vec=c(30, 50)){
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
  
  
  chosen_window = W_vec[1]
  
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
    }else{
      # on the previous steps data + new data, run a cvm test
      # How much data do we use? 
      
      if(num_active_batches < chosen_window){
        STAT <- cvm.res(x=data, y=f_ic)
        pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      }else{
        T_W = num_active_batches - chosen_window
        STAT <- cvm.res(x=data[T_W:num_active_batches], y=f_ic)
        pval <- find_cvm_pval_asym(n=m*chosen_window, STAT=STAT)
      }

      # find amount of previous data to remove
      if(pval < constant){
        chosen_window <- mean(c(chosen_window, W_vec[2]))
      }else{
        chosen_window = mean(c(chosen_window, W_vec[1]))
      }

      # num_remove <- prune(pval=pval, method=method)

      # for diagnostics
       num_keep <- chosen_window
      # store_keep <- append(store_keep, num_keep)
      # store_remove <- append(store_remove, num_remove)
      # print(paste("On Batch", i))
      # print(paste("Number removed equals ", num_remove))
      # print(paste("Amount of batches used equals ", num_keep))
      if(num_active_batches > W_vec[2]){
        data[,-1]
      }
    }
    # update control statistic 
    u_nK <- pval
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

Find_CL_CVM_Pval_DWindow_OOC <- function(ic_data=rnorm(500), 
                                      lambda=.05, 
                                      control_limit=.01, 
                                      m=5, 
                                      tau=0, 
                                      ICdist="rnorm", IC_dist_ops=NULL,
                                      OOCdist="rnorm", OOC_dist_ops=NULL,
                                      method="partha",
                                      W_vec=c(30, 50)){
  
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
  chosen_window = W_vec[1]
  
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
      
    }else{
      # on the previous steps data + new data, run a cvm test
      # How much data do we use? 
      
      if(num_active_batches < chosen_window){
        STAT <- cvm.res(x=data, y=f_ic)
        pval <- find_cvm_pval_asym(n=m*num_active_batches, STAT=STAT)
      }else{
        T_W = num_active_batches - chosen_window
        STAT <- cvm.res(x=data[T_W:num_active_batches], y=f_ic)
        pval <- find_cvm_pval_asym(n=m*chosen_window, STAT=STAT)
      }
      
      # find amount of previous data to remove
      if(pval < constant){
        chosen_window <- mean(c(chosen_window, W_vec[2]))
      }else{
        chosen_window = mean(c(chosen_window, W_vec[1]))
      }
      
      # num_remove <- prune(pval=pval, method=method)
      
      # for diagnostics
      # num_keep <- num_active_batches - num_remove
      # store_keep <- append(store_keep, num_keep)
      # store_remove <- append(store_remove, num_remove)
      # print(paste("On Batch", i))
      # print(paste("Number removed equals ", num_remove))
      # print(paste("Amount of batches used equals ", num_keep))
      if(num_active_batches > W_vec[2]){
        data[,-1]
      }
    }
    # update control statistic 
    u_nK <- pval
    u <- append(u, u_nK)
    i <- i + 1
    

  }
  i <- i - 1 # 10/20/2016 -- this adjust for the time out of control value!
  return(list("control_stats" = u, "Time OOC" = i))
}

