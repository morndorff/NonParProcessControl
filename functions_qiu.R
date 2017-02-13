# Qiu 2009 Functions

get_y <- function(x, boundaries, sum){
  int <- findInterval(x, boundaries)
  tab <- tabulate(int + 1, nbins=length(boundaries) + 1) # m*f
  emp_probs <- tab/length(x)  #f0
  return(emp_probs)
}


get_exact <- function(num_bps){
  exact_probs <- rep(1/num_bps, num_bps)
  return(exact_probs)
}

get_g <- function(data, quan, past_cols=0){
  #pastcols must be b/w 0 and dim(data)[2] - 1
  num_col <- dim(data)[2]
  num_quan <- length(quan)
  
  if(is.null(num_col) || num_col==1){
    if(past_cols!=0) stop("Too many columns on vector!")
    data_subset <- data
  }else{
    data_subset <- data[, (num_col - past_cols):num_col] 
  }
  int <- findInterval(data_subset, quan)
  tab <- tabulate(int + 1, nbins=num_quan + 1) # b/c needs positive, see ?tabulate for details
}

qiu_ARL <- function(ic_data=rnorm(500), kP=1, num_bps=5, control_limit, m=5, exact=FALSE, s=.01, ICdist="rnorm",
                    IC_dist_ops=NULL) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }
  
  if(exact){
    stop("Not Working Now")
    ic_probs <- get_exact(num_bps + 1)
    boundaries <- get_exact_boundaries(num_bps + 1)
    
  }else{
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0 # should be 1/(num_bps+1), approximately
  }

  data <- NULL
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  i <- 0
  num_bins <- num_bps + 1
  #print(mf0)
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops)) 
    # print("data")
    # print(data)
    g_n <- get_g(data = data, quan=boundaries, past_cols=0) #g(n)
    # print("boundaries")
    # print(boundaries)
    # print("gn")
    # print(g_n)
    g_n <- g_n + rnorm(num_bins, 0, s)
    # print(g_n)
    # C1 <- t((S_obs[, i] - S_exp[, i]) + (g_n -  mf0)) # check tranposition
    # C2 <- diag(1 / (S_exp[, i] + mf0))
    # C3 <- t(C1)
    # print(S_obs)
    # print(S_exp)
    # print(g_n)
    # print(mf0)
    C_n <- sum(((S_obs[, i] - S_exp[, i]) + (g_n -  mf0))^2 / (S_exp[, i] + mf0))
    
    
    # C_n <- C1 %*% C2 %*% C3
    # C_n <- as.vector(C_n)
    
    
    # print(paste("Cn", C_n))
    # print(paste("Ct", Ct))
    
    if(C_n <= kP){
      S_o_new = numeric(num_bps + 1)
      S_e_new = numeric(num_bps + 1)
      # print(paste("Cn < kP at time", i)) 
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kP) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kP) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    S_exp <- cbind(S_exp, S_e_new) 
    if(all(S_o_new==0)){
      u <- append(u, 0)
    }else{
      u_new <- sum((S_obs[, (i + 1)] - S_exp[, (i + 1)])^2 / (S_exp[, (i + 1)]))
      # print("U Test")
      # print(U_test)
      # U1 <- t(S_obs[, (i + 1)] - S_exp[, (i + 1)])
      # U2 <- diag( 1 / (S_exp[, (i + 1)]))
      # U3 <- t(U1)
      # print("Old U")
      # print(U1 %*% U2 %*% U3)
      u <- append(u, u_new)
    }
  }
  return(list("uP"=u, "Time OOC"=i))
}

qiu_Phase_II <- function(ic_data=rnorm(500), kP=1, num_bps=5, control_limit, m=5, exact=FALSE, 
                         tau=0, 
                         ICdist="rnorm", IC_dist_ops=NULL,
                         OOCdist="rnorm", OOC_dist_ops=NULL, s=.01){
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  # IC_dist_ops= NULL or list(mean=100, sd=2) or similar
  # ICdist/OOCdist= "qnorm" (for now, must be character
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }

  if(exact){
    stop("Not Working Now")
    ic_probs <- get_exact(num_bps + 1)
    boundaries <- get_exact_boundaries(num_bps + 1)
    
  }else{
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0 # should be 1/(num_bps+1), approximately
  }

  rIC <- get(ICdist, mode = "function", envir = parent.frame())

  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  data <- NULL

  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  i <- 0
  num_bins <- num_bps + 1

  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- do.call(rIC, c(list(m), IC_dist_ops))
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops))
    }
    g_n <- get_g(data = data, quan=boundaries, past_cols=0) #g(n)
    # print("boundaries")
    # print(boundaries)
    # print("gn")
    # print(g_n)
    g_n <- g_n + rnorm(num_bins, 0, s)
    # print(g_n)
    # C1 <- t((S_obs[, i] - S_exp[, i]) + (g_n -  mf0)) # check tranposition
    # C2 <- diag(1 / (S_exp[, i] + mf0))
    # C3 <- t(C1)
    # print(S_obs)
    # print(S_exp)
    # print(g_n)
    # print(mf0)
    C_n <- sum(((S_obs[, i] - S_exp[, i]) + (g_n -  mf0))^2 / (S_exp[, i] + mf0))
    
    
    # C_n <- C1 %*% C2 %*% C3
    # C_n <- as.vector(C_n)
    
    
    # print(paste("Cn", C_n))
    # print(paste("Ct", Ct))
    
    if(C_n <= kP){
      S_o_new = numeric(num_bps + 1)
      S_e_new = numeric(num_bps + 1)
      # print(paste("Cn < kP at time", i)) 
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kP) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kP) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    S_exp <- cbind(S_exp, S_e_new) 
    if(all(S_o_new==0)){
      u <- append(u, 0)
    }else{
      u_new <- sum((S_obs[, (i + 1)] - S_exp[, (i + 1)])^2 / (S_exp[, (i + 1)]))
      u <- append(u, u_new)
    }
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
}

qiu_L_ARL <- function(ic_data=rnorm(500), kL=1, num_bps=5, control_limit, m=5, exact=FALSE,
                      additive_constant=.001, s=.01) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }
  
  if(exact){
    ic_probs <- get_exact(num_bps + 1)
    boundaries <- get_exact_boundaries(num_bps + 1)
    
  }else{
    # Choose boundaries based on quantiles
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0
  }


  data <- NULL
  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  #print(mf0)
  i <- 0
  while(tail(u, 1) < control_limit) {
    i <- i + 1

    data <- cbind(data, rnorm(m))
    g_n <- get_g(data = data[, i], quan=boundaries, past_cols=0) #g(n)
    if(all(S_obs[, i]==0)){
      S_obs[, i] <- S_obs[, i] + additive_constant # avoids probs with log(0)
    }
    #print(mf0)
    C1 <- 2 * t(S_obs[, i] + g_n)
    C2 <- log( (S_obs[, i] + g_n) / (S_exp[, i] + mf0))
    #print(S_obs[, i])
    #print(S_exp[, i])
    #print(S_obs[, i] + g_n)
    #print(S_exp[, i] + mf0)
    #print(C1)
    #print(C2)
    C_n <- C1 %*% C2
    C_n <- as.vector(C_n)
    
    #print(C_n)
    if(C_n <= kL){
      S_o_new = numeric(num_bps + 1)
      S_e_new = numeric(num_bps + 1)
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kL) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kL) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    S_exp <- cbind(S_exp, S_e_new) 
    if(all(S_o_new==0)){
      u <- append(u, 0)
    }else{
      U1 <- t(S_obs[, (i + 1)])
      U2 <- log(S_obs[, (i+1)] / mf0)
      print(U1)
      print(U2)
      u_new <- 2 * U1 %*% U2
      u <- append(u, u_new)
    }
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
}

qiu_KS_ARL <- function(ic_data=rnorm(500), kK=.02, control_limit, m=5,
                       ICdist="rnorm", IC_dist_ops=NULL,
                       bootstrap_samples=3000) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input) (hK)
  # m: batch size
  # no exact option available here
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops)) 
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,  u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
  }
  return(list("uP"=u, "Time OOC"=i))
}


qiu_KS_PhaseII <- function(ic_data=rnorm(500), kK=.02, control_limit, m=5, exact=FALSE, 
                           tau=0, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL,
                           bootstrap_samples=3000,
                           keep_data=FALSE){
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){
    if(i <= tau){
      data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )

    }else{
      data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
    }
    f0_n <- fhat_ic(sort.int(data[, i]))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
    }else{
     if(i <= tau){
      data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 

    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK) 
    }
  }
  return(list("uP"=u, "Time OOC"=i))
}


qiu_CVM_ARL <- function(ic_data=rnorm(500), kK=.02, control_limit, m, bootstrap_samples=1000,
                        ICdist="rnorm", IC_dist_ops=NULL){
  
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- cvmtwo.res(ic_data, sample(ic_data, replace=T, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    
    D_n <- cvmtwo.res(data, ic_data)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
  
}

Estimate_Allowance <-
  function(ic_dist = "rnorm",
           M = 500,
           bootstrap_samples = 3000,
           tstat = SP_Two_Easy,
           m = 5,
           opt = "sd",
           runs = 6) {
    ic_dist <- get(ic_dist, mode = "function")
    if (opt == "sd")
      sd0 <- numeric(runs)
    if (opt == "mean")
      d0 <- numeric(runs)
    for (j in 1:runs) {
      ic_data <- ic_dist(M)
      D_n <- numeric(bootstrap_samples)
      for (i in 1:bootstrap_samples) {
        D <- tstat(ic_data, sample(ic_data, replace = T, size = m))
        D_n[i] <- D
      }
      if (opt == "sd") {
        sd0[j] <- sd(D_n)
      }
      if (opt == "mean") {
        d0[j] = mean(D_n)
      }
      
    }
    if (opt == "sd") {
      sd0 <- mean(sd0)
      #print("sd")
      return(sd0)
    }
    if (opt == "mean") {
      d0 <- mean(d0)
      #print("mean")
      return(d0)
    }
  }

Estimate_Allowance_Gamma <-
  function(ic_dist = "rnorm",
           M = 500,
           bootstrap_samples = 3000,
           tstat = cvm.res,
           m = 5,
           runs = 6) {
    ic_dist <- get(ic_dist, mode = "function")
    shape_est <- numeric(runs)
    rate_est <- numeric(runs)
    for (j in 1:runs) {
      ic_data <- ic_dist(M)
      D_n <- numeric(bootstrap_samples)
      for (i in 1:bootstrap_samples) {
        D <- tstat(ic_data, sample(ic_data, replace = T, size = m))
        D_n[i] <- D
      }
      params <- fitdistr(D_n, "gamma")
      shape_est[j] <- params$estimate["shape"]
      rate_est[j] <- params$estimate["rate"]
    }
    shape <- mean(shape_est)
    rate <- mean(rate_est)
    distribution_mean <- shape/rate
    return(c(shape, rate, distribution_mean))
  }






qiu_CUSUM_generic <- function(ic_data=rnorm(500), allowance_param=.02, control_limit, m, bootstrap_samples=1000,
                              ICdist="rnorm", IC_dist_ops=NULL, tstat){
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- tstat(x=ic_data, y=sample(ic_data, replace=F, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    
    D_n <- tstat(x=data, y=ic_data) # Use Test Statistic
    # this little x=, y= could cause problems later
    u_nK <- max(0,
                u[i] + (D_n - d0) - allowance_param)
    u <- append(u, u_nK)
  }
  return(list("uP"=u, "Time OOC"=i))
}
# Later, I'll need to add t-stat options, but for now....

# NOT ADDING TSTAT NULL CAUSE WE NEED IT
qiu_Phase_II_generic<- function(ic_data=rnorm(500), allowance_param=.02,  control_limit=20, m=5, exact=FALSE, 
                                tau=0, 
                                ICdist="rnorm", IC_dist_ops=NULL,
                                OOCdist="rnorm", OOC_dist_ops=NULL,
                                bootstrap_samples=2500, tstat){
  j <- 1:m
  # Not bootstrapping, will need to change language later
  # changed due to behavior of ties with SP_Test
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- tstat(ic_data, sample(ic_data, replace=F, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m 
  
  # Random variable generation functions
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- do.call(rIC, c(list(m), IC_dist_ops) ) 
      #print(data)
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    D_n <- tstat(data, ic_data)
    u_nK <- max(0,
                u[i] + (D_n - d0) - allowance_param)
    u <- append(u, u_nK)
  }
  return(list("uP"=u, "Time OOC"=i))
}


qiu_CVM_PhaseII <- function(ic_data=rnorm(500), allowance_parameter=.02,  control_limit=20, m=5, 
                            tau=0, 
                            ICdist="rnorm", IC_dist_ops=NULL,
                            OOCdist="rnorm", OOC_dist_ops=NULL,
                            bootstrap_samples=1000){
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- cvmtwo.res(ic_data, sample(ic_data, replace=T, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)

  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m 
  
  # Random variable generation functions
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- do.call(rIC, c(list(m), IC_dist_ops) ) 
      #print(data)
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    D_n <- cvmtwo.res(data, ic_data)
    u_nK <- max(0,
                u[i] + (D_n - d0) - allowance_parameter)
    u <- append(u, u_nK)
  }
  return(list("uP"=u, "Time OOC"=i))
}


R_EWMA_PhaseII <- function(lambda=.05, 
                           control_limit=.129375, 
                           m=5, 
                           tau=0, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL){
  
  
  
  # Initializing Variables  
  data <- NULL
  v <- 0 # Initialize v at 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(v, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- do.call(rIC, c(list(m), IC_dist_ops))
      #print(data)
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    
    x_bar <- mean(data)
    v_N <- lambda * x_bar + (1-lambda) * tail(v, 1)
    v <- append(v, v_N)
  }
  return(list("v_EWMA"=v, "Time OOC"=i))
}

R_EWMA_Find_CL <- function(lambda=.05, 
                           control_limit=.12, 
                           m=5, 
                           ICdist="rnorm", 
                           IC_dist_ops=NULL){
  # Initializing Variables  
  data <- NULL
  v <- 0 # Initialize v at 0
  i <- 0
  j <- 1:m
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(v, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    x_bar <- mean(data)
    v_N <- lambda * x_bar + (1-lambda) * tail(v, 1)
    v <- append(v, v_N)
  }
  return(list("v_EWMA"=v, "Time OOC"=i))
}

Find_ARL <- function(arl=200, lcl=0, ucl=.15, N_max=15, tol=1, sd_tol=3,
                          ICdist="rnorm", IC_dist_ops=NULL,
                          f=R_EWMA_Find_CL, N2_min=300, ...){
  # This can be adapted for a number of problems later
  # Works for any function where the output has "Time OOC"
  
  N <- 1
  current_lcl <- lcl
  current_ucl <- ucl
  arl_track <- arl + tol * 1.5 # will get overwritten later, jf 1st while iter
  
  while(N < N_max){
    
    # Is the arl found on the previous iteration near 200?
    # If yes, then output the list
    if(abs(arl_track - arl) < tol ){
      
      return(list("Calculated Control Limit"=new_cl, 
                  "Number of Iterations"=N,
                  "Calculated ARL"=arl_track,
                  "ARL_SD"=sd_arl))
    }
    
    
    new_cl <- mean(c(current_lcl, current_ucl)) # New control limit via bisection
    
    # Calculating control limit based on new control limit    
    # Two triggers: run for at least 500 iterations
    # and s.dev within sd_tol
    sd_arl <- sd_tol + 1
    arl_track <- NULL
    N2 <- 0
    while(sd_arl > sd_tol){
      new_arl <- f(control_limit=new_cl, ICdist = ICdist, IC_dist_ops = IC_dist_ops, ...)
      new_arl <- new_arl[["Time OOC"]]
      arl_track <- append(arl_track, new_arl)
      sd_arl <- sd(arl_track)/sqrt(length(arl_track))
      if(is.na(sd_arl)) sd_arl <- sd_tol + 1
      N2 <- N2 + 1
      if(N2 < N2_min) sd_arl <- sd_tol + 1 # don't stop until N2min
    }
    # output mean of arl_track
    arl_track <- mean(arl_track) # f(new_cl)
    print(paste("New Control Limit", new_cl, "has ARL of:", arl_track))
    print(paste("This took", N2, "iterations"))
    # Create new estimates
    if(arl_track < arl){
      current_lcl <- new_cl
    }else{
      current_ucl <- new_cl
    }
    
    N <- N + 1
  }
  stop("Did Not Converge")
}


EWMA_KS_Find_CL <- function(ic_data=rnorm(500),
                            lambda=.05, 
                            control_limit=.12, 
                            m=5, 
                            ICdist="rnorm", 
                            IC_dist_ops=NULL,
                            bootstrap_samples=3000){
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u<- 0 
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
    u <- append(u, u_nK)
    }
  return(list("u_EWMA_KS"=u, "Time OOC"=i))
}

EWMA_KS_PhaseII <- function(ic_data=rnorm(500),
                            lambda=.05, 
                           control_limit=0.033125, 
                           m=5, 
                           tau=0, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL,
                           bootstrap_samples=3000,
                           keep_data=FALSE){
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0 # Initialize u at 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){
      if(i <= tau){
        data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
      }else{
        data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
      }
      f0_n <- fhat_ic(sort.int(data[, i]))
      d_n1 <- f0_n - (j-1)/m 
      d_n2 <- (j/m) - f0_n
      D_n <- max(d_n1, d_n2)
      u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
      u <- append(u, u_nK)
    }else{
      if(i <= tau){
        data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        data <- do.call(rOOC, c(list(m), OOC_dist_ops))
      }
      f0_n <- fhat_ic(sort.int(data))
      d_n1 <- f0_n - (j-1)/m 
      d_n2 <- (j/m) - f0_n
      D_n <- max(d_n1, d_n2)
      u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
      u <- append(u, u_nK)
    }
  }
  return(list("u_EWMA"=u, "Time OOC"=i))
}
Find_ARL_MC <- function(arl=200, lcl=0, ucl=.15, N_max=15, tol=1, sd_tol=3,
                     ICdist="rnorm", IC_dist_ops=NULL,
                     f=R_EWMA_Find_CL, N2_min=300, ..., Cores=2){
  # This can be adapted for a number of problems later
  # Works for any function where the output has "Time OOC"
  #require(doMC)
  require(foreach)
  
  N <- 1
  current_lcl <- lcl
  current_ucl <- ucl
  arl_track <- arl + tol * 1.5 # will get overwritten later, jf 1st while iter
  registerDoMC(Cores)
  Num_Iter <- Cores * 5 #seems like a reasonable constant, don't want to add func arg
  while(N < N_max){
    
    # Is the arl found on the previous iteration near 200?
    # If yes, then output the list
    if(abs(arl_track - arl) < tol ){
      
      return(list("Calculated Control Limit"=new_cl, 
                  "Number of Iterations"=N,
                  "Calculated ARL"=arl_track,
                  "ARL_SD"=sd_arl))
    }
    
    
    new_cl <- mean(c(current_lcl, current_ucl)) # New control limit via bisection
    
    # Calculating control limit based on new control limit    
    # Two triggers: run for at least 500 iterations
    # and s.dev within sd_tol
    sd_arl <- sd_tol + 1
    arl_track <- NULL
    N2 <- 0
    while(sd_arl > sd_tol){
      par_obj <- foreach(n=1:Num_Iter) %dopar%{
        new_arl <- f(control_limit=new_cl, ICdist = ICdist, IC_dist_ops = IC_dist_ops, ...) 
        # if generic, make sure to include tstat in call
        new_arl <- new_arl[["Time OOC"]]
      }
      par_obj <- unlist(par_obj)
      arl_track <- append(arl_track, par_obj)
      sd_arl <- sd(arl_track)/sqrt(length(arl_track))
      if(is.na(sd_arl)) sd_arl <- sd_tol + 1
      N2 <- N2 + Num_Iter
      if(N2 < N2_min) sd_arl <- sd_tol + 1 # don't stop until N2min
    }
    print(paste("New Control Limit", new_cl, "has ARL of:", mean(arl_track)))
    print(paste("This took", length(arl_track), "iterations"))
    # output mean of arl_track
    arl_track <- mean(arl_track) # f(new_cl)

    # Create new estimates
    if(arl_track < arl){
      current_lcl <- new_cl
    }else{
      current_ucl <- new_cl
    }
    N <- N + 1
  }
  stop("Did Not Converge")
}

Find_ARL_MC_Reg <- function(arl=500, lcl=0, ucl=.15, N_max=15, sd_tol=3,
                        ICdist="rnorm", IC_dist_ops=NULL,
                        f=R_EWMA_Find_CL, N2_min=300, ..., Cores=2, num_points=5){
  # This can be adapted for a number of problems later
  # Works for any function where the output has "Time OOC"
  #require(doMC)
  require(foreach)
  
  registerDoMC(Cores)
  Num_Iter <- 5 * Cores
  x <- seq(lcl, ucl, length.out=num_points)
  y <- numeric(num_points)
  
  for(i in seq_along(x)){
    sd_arl <- sd_tol + 1
    arl_track <- NULL
    N2 <- 0
    while(sd_arl > sd_tol){
      par_obj <- foreach(n=1:Num_Iter) %dopar%{
        new_arl <- f(control_limit=x[i], ICdist = ICdist, IC_dist_ops = IC_dist_ops, ...)
        new_arl <- new_arl[["Time OOC"]]
      }
      par_obj <- unlist(par_obj)
      arl_track <- append(arl_track, par_obj)
      sd_arl <- sd(arl_track)/sqrt(length(arl_track))
      if(is.na(sd_arl)) sd_arl <- sd_tol + 1
      N2 <- N2 + Num_Iter
      if(N2 < N2_min) sd_arl <- sd_tol + 1 # don't stop until N2min
    }
    
    print(paste(" Control Limit", x[i], "has ARL of:", mean(arl_track)))
    print(paste("This took", length(arl_track), "iterations"))
    # output mean of arl_track
    y[i] <- mean(arl_track)
  }
  print(paste("X-Values are", x))
  print(paste("Y-Values are", y))
  
  m1 <- lm(x~y + I(y^2))
  print(arl)
  cl_pred <- predict(m1, newdata=data.frame(y=arl))
  cl_int <- predict(m1, newdata=data.frame(y=arl), interval="confidence")
  return(list("Prediction"=cl_pred, "Confidence Interval"=cl_int))
}



rmixnorm_1 <- function(N, u1, u2, sd1, sd2){
  # Does some mixture distributions
  components <- sample(1:2,prob=c(.5, .5),size=N,replace=TRUE)
  mus <- c(u1, u2)
  sds <- c(sd1, sd2)
  samples <- rnorm(n=N,mean=mus[components],sd=sds[components]) 
  return(samples)
}

rmixnorm <- function(N, u, s, probs) {
  if(length(u)!=length(s)) stop("Lengths don't match")
  len_param <- length(u)
  if(missing(probs)) probs <- rep(1/len_param, len_param)
  components <- sample(1:len_param, prob=probs, size=N, replace=TRUE)
  samples <- rnorm(n=N, mean=u[components], sd=s[components])
}


Ross_OOC_ARL = function(ic_data=rnorm(500), 
                        lambda=NULL, 
                        control_limit=NULL, 
                        m=5, 
                        tau=0, 
                        ICdist="rnorm", IC_dist_ops=NULL,
                        OOCdist="rnorm", OOC_dist_ops=NULL,
                        method="Ross",
                        ARL=200){
  require(cpm)
  if(method=="Ross") method="Cramer-von-Mises"
  M =length(ic_data)
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  i <- 1
  ic_tau = do.call(rIC, c(n=list(tau*m), IC_dist_ops))
  
  new_data = c(ic_data, ic_tau, do.call(rOOC, c(n=list(1000), OOC_dist_ops)))
  Detect_Time = detectChangePoint(new_data, cpmType="Cramer-von-Mises", ARL0=ARL, startup=M, lambda=NA) 
  Change_Detected=Detect_Time$changeDetected
  while(Change_Detected==FALSE){
    if(tau > 1000) stop("tau too big")
    new_data = c(new_data, do.call(rOOC, c(n=list(1000), OOC_dist_ops)))
    Detect_Time = detectChangePoint(new_data, cpmType=method, ARL0=ARL, startup=M, lambda=NA) 
    Change_Detected = Detect_Time$changeDetected
  }
  OOC_Time = ceiling((Detect_Time$detectionTime - M)/m) - tau
  return(list("Time OOC" = OOC_Time))
}


Find_N_CUSUM_ARL = function(ic_data = rnorm(500),
                            lambda = .1,
                            control_limit = 8,
                            m = 5,
                            ICdist = "rnorm",
                            IC_dist_ops = NULL,
                            bootstrap_samples = 3000,
                            track_candidates=NULL,
                            method="N-CUSUM"
                            ){
  if(method=="N-CUSUM"){
    ic_mean = mean(ic_data)
  }else if(method=="SD-CUSUM"){
    ic_sd = sd(ic_data)
  }
  ucl = control_limit
  lcl = -control_limit
  u_pos <- 0
  u_neg <- 0 # Initialize u at 0
  i <- 1
  
  while(tail(u_neg, 1) > lcl && tail(u_pos, 1) < ucl) {
    data <- do.call(ICdist, c(list(m), IC_dist_ops))
    
    if(method=="N-CUSUM"){
      data.mean = mean(data)
    
      u_pos_k <- max(0, tail(u_pos, 1) + (data.mean -ic_mean) - lambda)
      u_neg_k = min(0, tail(u_neg, 1) + (data.mean -ic_mean) + lambda)
    }else if(method=="SD-CUSUM"){
      data.sd = sd(data)
      
      u_pos_k <- max(0, tail(u_pos, 1) + (data.sd -ic_sd) - lambda)
      u_neg_k = min(0, tail(u_neg, 1) + (data.sd -ic_sd) + lambda)
      
    }
    u_pos <- append(u_pos, u_pos_k)
    u_neg <- append(u_neg, u_neg_k)
    i <- i + 1
  }
  # which triggered condition?

  i = i -1 #adjustment
  
  # If we have a large upper bound, we can calculate all our upper bounds below that number
  if(!is.null(track_candidates)){
    if(-ucl!=lcl) stop("this is not valid")
    track_ucl <-numeric(length(track_candidates))
    
    u = matrix(c(u_pos, u_neg) , ncol=length(u_pos), byrow=TRUE)
    u_abs = apply(u, 2, function(x) max(abs(x)))
    for(j in seq_along(track_candidates)){
      # for a particular candidate, see where things were triggered.
      # complication for upper and lower control limits
      #print(which.min(track_candidates[j] > u_pos))

      # track_upper[j] <- as.numeric(which.min(track_candidates[j] > u_pos)) - 1 #investigate later
      # if(track_upper[j] ==0) track_upper[j] = NA
      # track_lower[j] <- as.numeric(which.min(-track_candidates[j] < u_neg)) - 1 #investigate later
      # if(track_lower[j] ==0) track_lower[j] = NA
      
      track_ucl[j] = as.numeric(which.min(track_candidates[j] > u_abs)) - 1
      
      # track_ucl[j] = max(track_upper[j], track_lower[j], na.rm=TRUE)
      
      # track_ucl[j] <- as.numeric(which.min(track_candidates[j] > u)) - 1 #investigate later
    }

    names(track_ucl) <- as.character(track_candidates)
    return(list("Time OOC"=i, "Lower CLs"=track_ucl))
    
  }
  
  return(list("Time OOC" = i, "control_stats"= list("pos"=u_pos, "neg"=u_neg), "method"=method) )
}

CUSUM_OOC_ARL = function(ic_data=rnorm(500), 
                        lambda=NULL, 
                        control_limit=NULL, 
                        m=5, 
                        tau=0, 
                        ICdist="rnorm", IC_dist_ops=NULL,
                        OOCdist="rnorm", OOC_dist_ops=NULL,
                        method="N-CUSUM"){
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  if(method=="N-CUSUM"){
    ic_mean = mean(ic_data)
  }else if(method=="SD-CUSUM"){
    ic_sd = sd(ic_data)
  }
  
  ucl = control_limit
  lcl = -control_limit
  u_pos <- 0
  u_neg <- 0 # Initialize u at 0
  i <- 1

  while(tail(u_neg, 1) > lcl && tail(u_pos, 1) < ucl) {

    if(i <= tau){
      data <- do.call(rIC, c(list(m), IC_dist_ops))
    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    if(method=="N-CUSUM"){
      data.mean = mean(data)
      
      u_pos_k <- max(0, tail(u_pos, 1) + (data.mean -ic_mean) - lambda)
      u_neg_k = min(0, tail(u_neg, 1) + (data.mean -ic_mean) + lambda)
    }else if(method=="SD-CUSUM"){
      data.sd = sd(data)
      u_pos_k <- max(0, tail(u_pos, 1) + (data.sd -ic_sd) - lambda)
      u_neg_k = min(0, tail(u_neg, 1) + (data.sd -ic_sd) + lambda)
      
    }
    u_pos <- append(u_pos, u_pos_k)
    u_neg <- append(u_neg, u_neg_k)
    i <- i + 1
  }
  i = i -1 #adjustment
  
  return(list("Time OOC" = i, "control_stats"= list("pos"=u_pos, "neg"=u_neg)))
}



