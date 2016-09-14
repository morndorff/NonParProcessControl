# Effect of sample size on estimates

# Verdict: Doesn't matter much
newccfun <- function(num_candidates=30, 
                     lcl=.03, 
                     ucl=.055, 
                     bootstrap_samples=6000, 
                     ic_sample_size=500,
                     allowance_parameter=.05,
                     num_replications=10,
                     target_ARL=200,
                     runs=1000){
  ucl_mat <-matrix(, nrow=length(allowance_parameter), ncol=3)
  cl_mat <- matrix(nrow=2, ncol=length(lcl))
  cl_mat[1, ] = lcl
  cl_mat[2, ] = ucl
  rownames(cl_mat) <- c("lcl", "ucl")
  for(num in seq_along(allowance_parameter)){
    
    
    # Creating in control data
    ic_data <- rnorm(ic_sample_size)
    #ic_data <- rt(ic_sample_size, df=9)
    # Add noise to ic_data to get more accurate control limit
    
    # bootstrap sample from that in control data to get distribution of test statistic
    #in_control_tstat <- In_Control_Dist_CvM(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples)
    # in_control_tstat <- In_Control_Dist_CvM_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
    #                                             ICdist="rnorm", ic_params=list(mean=0, sd=1))
    in_control_tstat <- In_Control_Dist_CvM_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
                                                ICdist="rt", ic_params=list(df=5))
    
    # get mean of the distribution
    mean_ic_tstat <- mean(in_control_tstat)
    
    hist(in_control_tstat)
    plot(density(in_control_tstat))
    bw <- density(in_control_tstat)$bw
    
    
    # Upper run limit depends on what ARL you want to estimate
    upper_run_limit <- 5500 # go really high, but almost never fill in the matrix 
    control_limit_candidates <- seq(cl_mat["lcl",num], cl_mat["ucl", num], length.out=num_candidates)
    
    # Pre-allocating
    x <- numeric(upper_run_limit)
    cc <- matrix(nrow = runs, ncol=num_candidates)
    
    # Get control limit estimates by using the bootstrapped in control distribution
    # for the test statistic
    
    for(i in 1:runs){
      x <- numeric(upper_run_limit)
      j <- 2
      while(x[j-1] <= ucl[num]){
        means <- sample(in_control_tstat, 1)
        tstat_chosen <- rnorm(1, mean = means, sd = bw)
        x[j] <- allowance_parameter[num] * (tstat_chosen - mean_ic_tstat) + (1 - allowance_parameter[num])*x[j - 1]
        j <- j + 1
      }
      
      # for(j in 2:upper_run_limit){
      #   # For using kernel
      #   means <- sample(in_control_tstat, 1)
      #   tstat_chosen <- rnorm(1, mean = means, sd = bw)
      #   # tstat_chosen <- sample(in_control_tstat, 1)
      #   x[j] <- allowance_parameter * (tstat_chosen - mean_ic_tstat) + (1 - allowance_parameter)*x[j - 1]
      # }
      cc[i,] <- sapply(control_limit_candidates, 
                       function(control_limit_candidates) min(which(x > control_limit_candidates)))
      #print(i/runs)
    }
    # print(cc[1,])
    #print(test_means)
    cc_means <- colMeans(cc, na.rm=TRUE)
    #print(cc_means)
    names(cc_means) <- control_limit_candidates
    
    plot(control_limit_candidates, cc_means)
    print("done!")
    cc_means
    f <- approxfun(control_limit_candidates, cc_means)
    f_in <- approxfun(cc_means, control_limit_candidates)
    estimate_ucl <- f_in(200)
    ucl_mat[num, 1] <- estimate_ucl
    estimate_ucl <- f_in(180)
    ucl_mat[num, 2] <- estimate_ucl
    estimate_ucl <- f_in(250)
    ucl_mat[num, 3] <- estimate_ucl
  }
  colnames(ucl_mat) <- c("200", "180", "250")
  rownames(ucl_mat) <- as.character(allowance_parameter)
  ucl_mat
}

num_sims <- 20
num_samp_sizes <- 4
can_mat <- matrix(nrow=num_samp_sizes, ncol=num_sims)
samp_sizes <- c(100, 500, 5000, 10000)
for(i in 1:num_sims){
  can_mat[1,i] <- newccfun(ic_sample_size=samp_sizes[1])
  can_mat[2,i] <- newccfun(ic_sample_size=samp_sizes[2])
  can_mat[3,i] <- newccfun(ic_sample_size=samp_sizes[3])
  can_mat[4,i] <- newccfun(ic_sample_size=samp_sizes[4])
}

rowMeans(can_mat) # Conclusion: in control sample size has little effect
