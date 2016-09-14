# SPS Functions

Process_Stat <- function(proc, tau_minus, tstat, dist_ic, params, doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # tau_minus: A vector containing g(tau-,H0) (Should be of length tau-1)
  # tstat: A test statistic
  # dist_ic: The incontrol distribution - for wavelet and ks, should be "pnorm" or similar
  # params: parameters of the in control distribution
  lenproc <- dim(proc)[2] # run length of the process
  num.samp <- dim(proc)[1]
  # If the proc is initially a vector, its dim value will be 0, so we do this
  if(lenproc==1){
    ts <- NULL
    rdist_ic <- dist.conv.str(dist_ic, "r") # strip first letter off string and replace with "r"
    y <- get(rdist_ic, mode="function", envir=parent.frame())
    rvec <- do.call(y, c(list(num.samp), params))
    ts[1] <- do.call(tstat, c(list(rvec), list(dist_ic), params))
    ts[2] <- do.call(tstat, c(list(proc), list(dist_ic), params))
    ts[3] <- ts[2]-ts[1]
    STAT <- ts[3]
    tau_minus <- ts[1]
    which_max_tau <- 0
    return(list("Test Statistic"=STAT, "Tau Minus"=tau_minus, "Tau Estimate"= which_max_tau))
  }
  ts <- matrix(nrow=3,ncol=lenproc)
  # Handling g(tau-, H0)
  ts[1, 1:(lenproc-1)] <- tau_minus
  ts[1, lenproc] <- do.call(tstat, c(list(proc[,1:lenproc]), list(dist_ic), params))
  
  if(lenproc!=1) {
  for(tau in 0:(lenproc-1)){
    ts[2, (tau+1)] <- do.call(tstat, c(list(proc[, (tau+1):lenproc]), list(dist_ic), params)) #g(H0, tau+)
  }
  ts[3, ] <- ts[2, ] - ts[1, ] # Take difference
  }
  STAT <- max(ts[3,])
  which_max_tau <- which.max(ts[3,])-1
  tau_minus <- ts[1,]
  if(detail){
    deta <- list("Test Statistic"=STAT,"All Stats"=ts, "Tau Minus"=tau_minus, "Tau Estimate"=which_max_tau)
    return(deta)
  }
  return(list("Test Statistic"=STAT, "Tau Minus"=tau_minus, "Tau Estimate"=which_max_tau))
}

Find_IC_RL_Slow <- function(num.samp=32, 
                            dist="pnorm", params=list(mean=0, sd=1), 
                            tstat=wave.energy, UCL){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  count <- 1
  track_stat <- NULL
  tau_minus <- NULL
  tau_est <- NULL
  tstat_iter <- 0
  tstat_proc <- 0
  dist_ic <- dist
  dist <- substring(dist,first=2)
  while(tstat_iter < max(UCL)){
    Proc <- Sim_IC_Process_Iter(proc=Proc, num.samp=num.samp, dist=dist, 
                                param=params) # good
    tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                               dist_ic=dist_ic, params=params, tau_minus=tau_minus)
    tau_minus <- tstat_proc[["Tau Minus"]]
    tstat_iter <- tstat_proc[["Test Statistic"]]
    track_stat <- append(track_stat, tstat_iter)
    count <- count +1
  }
  tau_est <- tstat_proc[["Tau Estimate"]]
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "RL for Corresponding UCL"=RL, 
              "UCLS"=UCL, "Estimated Tau"=tau_est)
  return(res)
}

Find_CP_RL_Slow <- function(num.samp=32, dist_one="pnorm", param_one=list(mean=0, sd=1),
                            dist_two="pnorm", param_two=list(mean=0, sd=2), cp=1,
                            tstat=wave.energy, UCL){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  tau_minus <- NULL
  track_stat <- NULL
  tstat_proc <- 0
  tstat_iter <- 0
  dist_ic <- dist_one
  dist_one <- substring(dist_one, first=2)
  dist_two <- substring(dist_two, first=2)  
  while(tstat_iter < max(UCL)){
    Proc <- Sim_CP_Process_Iter(proc=Proc, num.samp=num.samp, cp=cp, 
                                dist_one=dist_one, param_one=param_one,
                                dist_two=dist_two, param_two=param_two)# good
    tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                               dist_ic=dist_ic, params=param_one, tau_minus=tau_minus)
    tau_minus <- tstat_proc[["Tau Minus"]]
    tstat_iter <- tstat_proc[["Test Statistic"]]
    track_stat <- append(track_stat, tstat_iter)
  }
  tau_est <- tstat_proc[["Tau Estimate"]]
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "RL for Corresponding UCL"=RL, 
              "UCLS"=UCL, "Estimated Tau"=tau_est)
  return(res)
}

ARL_Proc_Tol <- function(UCL, num.samp = 30,
                          tolerance=10,
                          method=Find_IC_RL_Fast, 
                          tstat=wave.energy,
                          dist="pnorm",
                          params=list(mean=0, sd=1),
                          min_rl=10,
                          ...){
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  tol <- tolerance +1
  e_time <- 0
  count <- 0
  len_UCL <- length(UCL)
  track_max_rl <- NULL
  # Calculating ARLs
  while(tol > tolerance){
    RLs_det <- method(num.samp=num.samp, 
                      dist=dist, 
                      params=params,
                      tstat=tstat, 
                      UCL=UCL, ...)
    RLs[[length(RLs)+1]] <- RLs_det[["RL for Corresponding UCL"]]
    Max_RL <- max(unlist(RLs_det[["RL for Corresponding UCL"]]))  
    track_max_rl <- append(track_max_rl,Max_RL)
    if(length(track_max_rl)>min_rl)  tol <- sd(track_max_rl)/sqrt(length(track_max_rl))
    howlong <- proc.time()-ptm
    e_time <- howlong["elapsed"]
    count <- count +1
    print(paste("Iteration ",count, ", Total Time: ",e_time, " RL:", Max_RL))
  }
  
  matRL <- matrix(unlist(RLs), ncol=len_UCL, byrow=TRUE)# Making ARL Matrix
  ARL <- matrix(,nrow=2, ncol=len_UCL) 
  ARL[1,] <- colMeans(matRL)
  ARL[2,] <- apply(matRL, 2, sd)
  colnames(ARL) <- as.character(round(UCL,3))
  rownames(ARL) <- c("mean", "sd")
  print(ARL)
  return(list(ARL, RLs, e_time))
}


ARL_Proc <- function(UCL, num.samp = 30,
                     time=60, 
                     method=Find_IC_RL_Fast, 
                     tstat=wave.energy,
                     dist="pnorm",
                     params=list(mean=0, sd=1), 
                     ...){
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  e_time <- 0
  count <- 0
  len_UCL <- length(UCL)
  # Calculating ARLs
  while(e_time < time){
    RLs_det <- method(num.samp=num.samp, 
                      dist=dist, 
                      params=params,
                      tstat=tstat, 
                      UCL=UCL, ...)
    RLs[[length(RLs)+1]] <- RLs_det[["RL for Corresponding UCL"]]
    Max_RL <- max(unlist(RLs_det[["RL for Corresponding UCL"]]))    
    howlong <- proc.time()-ptm
    e_time <- howlong["elapsed"]
    count <- count +1
    print(paste("Iteration ",count, ", Total Time: ",e_time, " RL:", Max_RL))
  }
  
  matRL <- matrix(unlist(RLs), ncol=len_UCL, byrow=TRUE)# Making ARL Matrix
  ARL <- matrix(,nrow=2, ncol=len_UCL) 
  ARL[1,] <- colMeans(matRL)
  ARL[2,] <- apply(matRL, 2, sd)
  colnames(ARL) <- as.character(round(UCL,3))
  rownames(ARL) <- c("mean", "sd")
  print(ARL)
  return(list(ARL, RLs, e_time))
}

Find_ARL_OOC <- function(num.samp=30, dist_one="qnorm", param_one=list(mean=0, sd=1),
                         dist_two="qnorm", param_two=list(mean=0, sd=2), cp=1,
                         tstat=wave.energy, UCL, time = 30, method=Find_CP_RL_Fast, ...){
  # Dist_one and dist_two need to be either qnorm or pnorm
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  Tau_Estimate <- vector(mode="list", length=0)
  e_time <- 0
  len_UCL <- length(UCL)
  count <- 0
  # Calculating ARLs
  while(e_time < time){
    RLs_det <- method(num.samp=num.samp, 
                      dist_one=dist_one, param_one=param_one,
                      dist_two=dist_two, param_two=param_two,
                      tstat=tstat, UCL=UCL, cp=cp, ...)
    RLs[[length(RLs)+1]] <- RLs_det[["RL for Corresponding UCL"]] # Append list of RL's
    Tau_Estimate[[length(Tau_Estimate)+1]]<- RLs_det[["Estimated Tau"]]
    howlong <- proc.time()-ptm
    e_time <- howlong["elapsed"]
    count <- count +1
    print(paste("Iteration ",count, ", Total Time: ",e_time))
  }
  Tau_Mean <- mean(unlist(Tau_Estimate))
  Tau_SD <- sd(unlist(Tau_Estimate))#/sqrt(length(unlist(Tau_Estimate)))
  matRL <- matrix(unlist(RLs), ncol=len_UCL, byrow=TRUE)# Making ARL Matrix
  ARL <- matrix(,nrow=2, ncol=len_UCL) 
  ARL[1,] <- colMeans(matRL)
  ARL[2,] <- apply(matRL, 2, sd)
  colnames(ARL) <- as.character(round(UCL,3))
  rownames(ARL) <- c("mean", "sd")
  print(ARL)
  return(list(ARL, RLs, e_time, "Tau Estimate"=Tau_Estimate, "Avg Tau"=Tau_Mean, "Tau SD"=Tau_SD))
}

Process_Stat_W <- function(proc, tau_minus, tstat, 
                           WSize, dist_ic, params, 
                           throw=TRUE,doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # tau_minus: A vector containing g(tau-,H0) (Should be of length tau-1)
  # tstat: A test statistic
  # dist_ic: The incontrol distribution - for wavelet and ks, should be "pnorm" or similar
  # params: parameters of the in control distribution
  lenproc <- dim(proc)[2] # run length of the process
  num.samp <- dim(proc)[1]
  # Should never be needed, unless Window Size is 1
  if(lenproc==1){
    if(WSize > 1) stop("This should never occur unless WSize=1")
    ts <- NULL
    # ts[1] <- 0
    rdist_ic <- dist.conv.str(dist_ic, "r") # Convert from pnorm to rnorm
    y <- get(rdist_ic, mode="function", envir=parent.frame())
    rvec <- do.call(y, c(list(num.samp), params))
    ts[1] <- do.call(tstat, c(list(rvec), list(dist_ic), params))
    ts[2] <- do.call(tstat, c(list(proc), list(dist_ic), params))
    #ts[2] <- tstat(proc, dist_ic, ...)
    ts[3] <- ts[2]-ts[1]
    STAT <- ts[3]
    tau_minus <- ts[1]
    which_max_tau <- 0
    return(list("Test Statistic"=STAT, "Tau Minus"=tau_minus, "Tau Estimate"=which_max_tau))
  }
  
  ts <- matrix(nrow=3,ncol=WSize) 
  # NOTE: under option 1, we have to recalculate tau- at each point
  # NOTE: under option 2, we do NOT
  # Proceeding under option 1
  
  # Handling g(tau-, H0)
  # ts[1, 1:(lenproc-1)] <- tau_minus
  # ts[1, lenproc] <- do.call(tstat, c(list(proc[,1:lenproc]), list(dist_ic), params))
  if(lenproc!=1) {
    T_minus_W <- lenproc-WSize
    tr <-1 
    for(tau in T_minus_W:(lenproc-1)){
      if(throw){
        ts[1, tr] <- do.call(tstat, c(list(proc[, T_minus_W:tau]), list(dist_ic), params)) #g(H0, tau-)
        
      } else {
        ts[1, tr] <- do.call(tstat, c(list(proc[, 1:tau]), list(dist_ic), params)) #g(H0, tau-)
      }
      ts[2, tr] <- do.call(tstat, c(list(proc[, (tau+1):lenproc]), list(dist_ic), params)) #g(H0, tau+)
      #print(tau)
      tr <- tr+1
    }
    ts[3, ] <- ts[2, ] - ts[1, ] # Take difference
  }
  #print(ts)
  STAT <- max(ts[3,])
  which_max_tau <- T_minus_W + which.max(ts[3,]) -1
  if(detail){
    deta <- list("Test Statistic"=STAT, "All Stats"=ts, "Tau Estimate"=which_max_tau)
    return(deta)
  }
  return(list("Test Statistic"=STAT, "Tau Estimate"=which_max_tau))
}

Find_IC_RL_Windowed <- function(num.samp=30, 
                                dist="pnorm", params=list(mean=0, sd=1), 
                                tstat=wave.energy, UCL, WSize=10, throw=TRUE, doplot=FALSE){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  count <- 1
  track_stat <- NULL
  tstat_proc <- 0
  tstat_iter <- 0
  dist_ic <- dist
  dist <- substring(dist,first=2)
  tau_minus <- NULL
  while(tstat_iter < max(UCL)){
    Proc <- Sim_IC_Process_Iter(proc=Proc, num.samp=num.samp, dist=dist, 
                                param=params) # good
    lenproc <- dim(Proc)[2]
    # If W >= to T, call original function
    if(WSize >= lenproc){
      tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                                 dist_ic=dist_ic, params=params, tau_minus=tau_minus)
      tau_minus <- tstat_proc[["Tau Minus"]]
      tstat_iter <- tstat_proc[["Test Statistic"]]
      track_stat <- append(track_stat, tstat_iter)
      count <- count +1
    } else { # If W > T, calculate using Windowed Version
      tstat_proc <- Process_Stat_W(proc=Proc, tstat=tstat, 
                                   dist_ic=dist_ic, params=params, WSize=WSize, throw=throw)
      tstat_iter <- tstat_proc[["Test Statistic"]]
      track_stat <- append(track_stat, tstat_iter)
      #track_stat <- append(track_stat, tstat_proc)
      count <- count +1
    }
  }
  tau_est <- tstat_proc[["Tau Estimate"]]
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  if(doplot) plot(track_stat)
  res <- list("h(t)"=track_stat, "RL for Corresponding UCL"=RL, "UCLS"=UCL, "Estimated Tau"=tau_est)
  return(res)
}

Find_CP_RL_Windowed<- function(num.samp=30, dist_one="pnorm", param_one=list(mean=0, sd=1),
                               dist_two="pnorm", param_two=list(mean=0, sd=2), cp=10,
                               tstat=wave.energy, UCL, WSize=10, throw=TRUE){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  track_stat <- NULL
  tstat_proc <- 0
  tstat_iter <- 0
  tau_minus <- NULL
  dist_ic <- dist_one
  dist_one <- substring(dist_one, first=2)
  dist_two <- substring(dist_two, first=2)  
  while(tstat_iter  < max(UCL)){
    Proc <- Sim_CP_Process_Iter(proc=Proc, num.samp=num.samp, cp=cp, 
                                dist_one=dist_one, param_one=param_one,
                                dist_two=dist_two, param_two=param_two)# good
    lenproc <- dim(Proc)[2]
    # If W >= to T, call original function
    if(WSize >= lenproc){
      tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                                 dist_ic=dist_ic, params=param_one, tau_minus=tau_minus)
      tau_minus <- tstat_proc[["Tau Minus"]]
      tstat_iter <- tstat_proc[["Test Statistic"]]
      track_stat <- append(track_stat, tstat_iter)
    } else { # If W > T, calculate using Windowed Version
      tstat_proc <- Process_Stat_W(proc=Proc, tstat=tstat, 
                                   dist_ic=dist_ic, params=param_one, 
                                   WSize=WSize, throw=throw)
      #print(paste("Used Window, time=", lenproc))
      tstat_iter <- tstat_proc[["Test Statistic"]]
      track_stat <- append(track_stat, tstat_iter)
    }
  }
  tau_est <- tstat_proc[["Tau Estimate"]]
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "RL for Corresponding UCL"=RL,
              "UCLS"=UCL, "Estimated Tau"=tau_est)
  return(res)
}

Sim_IC_Process <- function(num.samp, run.length, dist, params){
  # Simulates a process of run length N with samples of size k
  # Inputs:
  # num.samp: The size of the samples. 1 corresponds to one data point per run unit
  # run.length: Length of the run
  # dist: what distribution do the data come from
  # params: what parameters to call from the dist function
  total_samps <- num.samp*run.length
  Proc <- matrix(make_sample(total_samps,dist,params),nrow=num.samp,ncol=run.length)
}

Sim_IC_Process_Iter <- function(proc=NULL, num.samp, dist, params){
  # Simulate a Incontrol Process Iteratively
  if(is.null(proc)){
    proc <- matrix(make_sample(num.samp, dist, params), nrow=num.samp, ncol=1)
    return(proc)
  }
  new_samp <- make_sample(num.samp, dist, params)
  proc <- cbind(proc, new_samp)
}

Sim_CP_Process_Iter <- function(proc=NULL, num.samp, cp, dist_one, param_one,
                                dist_two, param_two){
  # Args:
  # dist_one: "norm" "unif"
  # param_one: list(mean=0, sd=1)
  # cp: When does the distribution change.  If it is one, it does so immediately
  if(is.null(proc)){
    if(cp==1){
      proc <- matrix(make_sample(n=num.samp, dist=dist_two, params=param_two), nrow=num.samp, ncol=1)
      return(proc)     
    }
    proc <- matrix(make_sample(n=num.samp, dist=dist_one, params=param_one), nrow=num.samp, ncol=1)
    return(proc)
  }
  lenproc <- dim(proc)[2]
  if(lenproc < cp){
    new_samp <- make_sample(n=num.samp, dist=dist_one, params=param_one)
    proc <- cbind(proc, new_samp)
    return(proc)
  }
  if(lenproc >= cp){
    new_samp <- make_sample(n=num.samp, dist=dist_two, params=param_two)
    proc <- cbind(proc, new_samp)
    return(proc)
  }
}

Sim_CP_Process <- function(num.samp,run.length,dist_one, param_one, 
                           dist_two, param_two, bpoint){
  samp_one <- (num.samp*bpoint)
  samp_two <- (num.samp*(run.length-bpoint))
  Proc_One <- matrix(make_sample(samp_one, dist_one, param_one),
                     nrow=num.samp, ncol=bpoint)
  Proc_Two <- matrix(make_sample(samp_two, dist_two, param_two),
                     nrow=num.samp, ncol=run.length-bpoint)
  CP_Proc <- cbind(Proc_One, Proc_Two)
  attr(CP_Proc, "bp") <- bpoint
  return(CP_Proc)
}



Find_IC_RL_Fast <- function(num.samp=32, 
                            dist="pnorm", params=list(mean=0, sd=1), 
                            tstat=wave.energy,
                            UCL, detail=FALSE, weight=FALSE){
  # num.samp - numeric, # size of RS from each dist
  # dist - Incontrol Distribution - prefer string
  # params - Incontrol Distribution Parameters -- MUST BE LIST
  # tstat - Test Statistic to use, i.e. ks.res.simp, wave.den
  # UCL - vector of potential upper control limits -- max controls ends of loop
  # detail - Gives more output
  # weight - Weighting scheme used. Currently if weight=FALSE, then uses sample mean
  Proc <- NULL
  h_t <- 0
  h_t_new <- 0
  g_t_p <- NULL
  g_t_m <- NULL
  which_tau <- NULL
  which_tau_weight <- NULL
  dist_ic <- dist
  dist <- substring(dist,first=2)
  while(h_t_new < max(UCL)){
    Proc <- Sim_IC_Process_Iter(proc=Proc, num.samp=num.samp, dist=dist, 
                                param=params) # New Realization @ Time T
    Time <- dim(Proc)[2] # Can probably just be replaced by count var
    g_t <- update_g_t(nvec=Proc[, Time], 
                      theta_p = g_t_p, 
                      theta_m = g_t_m, 
                      tstat = tstat, 
                      dist_ic = dist_ic,
                      params = params) # get new estimates for g+ and g-
    g_t_p <- g_t[, 1]
    g_t_m <- g_t[, 2]
    # Scaling
    if(weight==TRUE){
      pos_scale <- seq(1/Time, 1, length.out=Time)
      neg_scale <- rev(pos_scale)
    }
    if(weight==FALSE){
      pos_scale <- rev(1 / 1:Time)
      neg_scale <- 1 / 1:Time
    }
    h_g_p <- g_t_p * pos_scale 
    h_g_m <- g_t_m * neg_scale
    which_tau <- append(which_tau, which.max(h_g_p - h_g_m)) # Diagnostics for Tau selection
    which_tau_weight <- append(which_tau_weight, tail(which_tau/Time,1))
    h_t <- append(h_t, max(h_g_p - h_g_m)) # For finding RLs
    h_t_new <- tail(h_t, 1) # When does loop stop?
  }
  RL <- sapply(UCL, function(x) min(which(h_t >= x))) # Finding RLs to Corresponding UCL
  lenproc <- dim(Proc)[2]
  if(detail){
    res <- list("h_t"=h_t, 
                "Length of Process"=lenproc, 
                "RL for Corresponding UCL"=RL, 
                "UCLs"=UCL,
                "hgp"=h_g_p,
                "hmp"=h_g_m,
                "Tau Chosen at each Time Period"=which_tau,
                "Tau Chosen Divided by Time" = which_tau_weight)
    return(res)
  }
  res <- list("h_t"=h_t, "Length of Process"=lenproc, "RL for Corresponding UCL"=RL, "UCLs"=UCL)
  return(res)
}

Find_CP_RL_Fast <- function(num.samp=32, dist_one="pnorm", param_one=list(mean=0, sd=1),
                            dist_two="pnorm", param_two=list(mean=0, sd=2), cp=1,
                            tstat=wave.energy, UCL, detail=FALSE, weight=TRUE){
  Proc <- NULL
  h_t <- 0
  h_t_new <- 0
  g_t_p <- NULL
  g_t_m <- NULL
  which_tau <- NULL
  which_tau_weight <- NULL
  dist_ic <- dist_one
  dist_one <- substring(dist_one, first=2)
  dist_two <- substring(dist_two, first=2)  
  while(h_t_new < max(UCL)){
    Proc <- Sim_CP_Process_Iter(proc=Proc, num.samp=num.samp, cp=cp, 
                                dist_one=dist_one, param_one=param_one,
                                dist_two=dist_two, param_two=param_two) # New Realization @ Time T
    Time <- dim(Proc)[2] # Can probably just be replaced by count var
    g_t <- update_g_t(nvec=Proc[, Time], 
                      theta_p = g_t_p, 
                      theta_m = g_t_m, 
                      tstat = tstat, 
                      dist_ic = dist_ic,
                      params=param_one) # get new estimates for g+ and g-
    g_t_p <- g_t[, 1]
    g_t_m <- g_t[, 2]
    if(weight==TRUE){
      pos_scale <- seq(1/Time, 1, length.out=Time)
      neg_scale <- rev(pos_scale)
    }
    if(weight==FALSE){
      pos_scale <- rev(1 / 1:Time)
      neg_scale <- 1 / 1:Time
    }
    h_g_p <- g_t_p * pos_scale
    h_g_m <- g_t_m * neg_scale
    which_tau <- append(which_tau, which.max(h_g_p - h_g_m))
    which_tau_weight <- append(which_tau_weight, tail(which_tau/Time,1))
    h_t <- append(h_t, max(h_g_p - h_g_m))
    #h_t <- append(h_t, find_max_dif(theta)) # New estimates for h(t)
    h_t_new <- tail(h_t, 1) # For finding RLs
  }
  RL <- sapply(UCL, function(x) min(which(h_t >= x))) # Finding RLs to Corresponding UCL
  lenproc <- dim(Proc)[2]
  if(detail){
    res <- list("h_t"=h_t, 
                "Length of Process"=lenproc, 
                "RL for Corresponding UCL"=RL, 
                "UCLs"=UCL,
                "hgp"=h_g_p,
                "hmp"=h_g_m,
                "Tau Chosen at each Time Period"=which_tau,
                "Tau Chosen Divided by Time" = which_tau_weight)
    return(res)
  }
  res <- list("h_t"=h_t, "Length of Process"=lenproc, "RL for Corresponding UCL"=RL, "UCLs"=UCL)
  return(res)
}

update_g_t <- function(nvec, theta_p, theta_m, tstat, dist_ic, params, ...){
  # Updates 'Fast' Algorithm for changepoint process
  # Input:
  # rlength: T, time
  # nvec: new vector recieved at time T
  # theta_p tau+ estimates
  # theta_m tau- estimates
  # output: 2 column matrix of updated parameters 
  rlength <- length(theta_p) +1 
  theta_p <- append(theta_p, 0)
  theta_m <- append(theta_m, 0)# 
  nvec_null <- do.call(tstat, c(list(nvec), list(dist_ic), params))
  #nvec_null <- tstat(nvec, dist_ic, ...)
  if(rlength==1){
    # Special case when T=1
    # TODO: Averaging
    # Generate random sample from null distribution here, test against null
    rdist_ic <- dist.conv.str(dist_ic, "r") # Convert from pnorm/qnorm to rnorm
    y <- get(rdist_ic, mode="function", envir=parent.frame())
    rvec <- y(length(nvec), ...)
    theta_m <- tstat(rvec, dist_ic, ...)
    theta_p <- nvec_null
    return(cbind(theta_p,theta_m))
  }
  theta_m[rlength] <- theta_m[(rlength-1)] + theta_p[(rlength-1)] # Will Break if T=1
  theta_p <- theta_p + nvec_null
  return(cbind(theta_p,theta_m))
}


find_max_dif <- function(theta){
  # NOT Absolute Values
  h_t <- max(theta[,1]-theta[,2])
  return(h_t)
}