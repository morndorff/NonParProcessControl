# EWMA vs CUSUM for control charts
#CL_Finding_With_UB(allow_param=.05, cl_seq=c(.055, .05, .048, .045), f_in_control=EWMA_Find_CL_CVM, sd_tol_init=30, Cores=4)
# Still needs a little work, but this works well for now

rm(list=ls())
setwd("/home/mark/Dropbox/Research/GoF Test/")
source("functions.R")
set.seed(5000)
library(doMC)
library(foreach)

# Params 
M <- 500
rdist <- rnorm
rdist_character <- "rnorm"
ic_params <- list(mean=0, sd=1)

# Adjusted Distribution
adist <- rt
adist_character <- "rt"
a_params = list(df=5)

ARL <- 200
margin_error <- 13
m <- 5
Cores <- 4
registerDoMC(Cores)

# Find Control Control Limits for the Process
##### Parameters for UCL/LCL
# Most Important
f_in_control <- EWMA_Find_CL_CVM
f_out_control <- EWMA_CVM_OOC
tstat <- cvm.res # An available one sample test statistic
lambdas = c(.05, .1, .2, .3)



# Writing new function for multiple allowance parameters

newccfun <- function(num_candidates=30, 
                     lcl=.03, 
                     ucl=.055, 
                     bootstrap_samples=6000, 
                     ic_sample_size=500,
                     allowance_parameter=.05,
                     num_replications=10,
                     target_ARL=200,
                     runs=1000,
                     df_opt=TRUE){
  ucl_mat <-matrix(, nrow=length(allowance_parameter), ncol=3)
  cl_mat <- matrix(nrow=2, ncol=length(lcl))
  cl_mat[1, ] = lcl
  cl_mat[2, ] = ucl
  rownames(cl_mat) <- c("lcl", "ucl")
  for(num in seq_along(allowance_parameter)){
    estimate_ucl_200 <- numeric(num_replications)
    estimate_ucl_170 <- numeric(num_replications)
    estimate_ucl_250 <- numeric(num_replications)
    
    for(reps in 1:num_replications){
      
      # Creating in control data
      ic_data <- rnorm(ic_sample_size)
      #ic_data <- rt(ic_sample_size, df=9)
      # Add noise to ic_data to get more accurate control limit
      
      # bootstrap sample from that in control data to get distribution of test statistic
      #in_control_tstat <- In_Control_Dist_CvM(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples)
      # in_control_tstat <- In_Control_Dist_CvM_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
      #                                             ICdist="rnorm", ic_params=list(mean=0, sd=1))
      if(df_opt){
        in_control_tstat <- In_Control_Dist_CvM_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
                                                    ICdist="rt", ic_params=list(df=5*allowance_parameter[num] / .07)) ## HARD CODE CITY, CAPTAIN 
        
      }else{
      in_control_tstat <- In_Control_Dist_CvM_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
                                                  ICdist="rt", ic_params=list(df=5))
      }
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
      estimate_ucl_200[reps] <- f_in(200)
      print(estimate_ucl_200)
      # ucl_mat[num, 1] <- estimate_ucl_200
      estimate_ucl_170[reps] <- f_in(170)
      #ucl_mat[num, 2] <- estimate_ucl_180
      estimate_ucl_250[reps] <- f_in(250)
      #ucl_mat[num, 3] <- estimate_ucl_250
    }
    ucl_mat[num, 1] <- mean(estimate_ucl_200)
    ucl_mat[num, 2] <- mean(estimate_ucl_170)
    ucl_mat[num, 3] <- mean(estimate_ucl_250)
    
  }
  colnames(ucl_mat) <- c("200", "170", "250")
  rownames(ucl_mat) <- as.character(allowance_parameter)
  ucl_mat
}





control_limit_mat <- newccfun(ic_sample_size=M, allowance_parameter = lambdas, lcl=c(.03, .07, .12, .17), ucl=c(.06, .11, .185, .25))
control_limit_mat
is.na(control_limit_mat)


# Get initial approximations
# mat_init_cl <- Initial_CL_Approx_CVM(num_candidates = 30,
#                       ARL = ARL,
#                       bootstrap_samples = 5000,
#                       ic_sample_size = M,
#                       allowance_parameter = lambdas,
#                       ic_dist=rdist_character,
#                       ic_params=ic_params,
#                       runs = 1000,
#                       avg_runs = 10,
#                       est_high=TRUE,
#                       over_estimation_param=1.6,
#                       m=m)
lcl <- control_limit_mat[, "170"]
ucl <- control_limit_mat[, "250"]


num_lambda <- length(lambdas)

# Tolerance Parameters
sd_tol <- 10
N2 <- 0 
N2_min <- 60 # At least 300 runs
N2_max <- 10000
sd_arl <- sd_tol + 2

arl_track <- numeric(0)
sd_track <- numeric(0)

cl_pred <- numeric(num_lambda)
cl_int <- matrix(nrow=num_lambda, ncol=3) #lower and upper confidence band
colnames(cl_int) <- c("fit", "lwr", "upr")

num_points <- 5 # For regression estimate
x <- matrix(nrow = num_lambda, ncol = num_points)
y <- matrix(nrow = num_lambda, ncol = num_points)
for(i in 1: num_lambda){
  x[i , ] <- seq(lcl[i], ucl[i], length.out=num_points)
}
# Fill in x matrix
temp_results <- vector(mode="list", length=length(lambdas))
for(i in 1:num_lambda){
  temp_results[[i]] <-
    Iterative_CL_Finding_MC(allow_param = lambdas[i],
                         cl_seq=x[i, ],
                         f_in_control = f_in_control,
                         margin_error = margin_error)
  saveresults(paste0("Iteration", i, ".RData"))
}


# Predictions for Control Limits

cl_pred <- sapply(temp_results, function(x) x[["Predicted CL"]])
saveresults("halfway_done")


allow_params <- lambdas
num_allow_params <- length(lambdas)

mean_shifts <- c(1, .9, .8, .7 , .6, .5, .4, .3, .2)
sd_tol <- c(.25, .3, .4, .45, .5, 2, 4, 6, 8)  # temporary
ooc_arl <- numeric(0)
ARLs <- list(num_allow_params)
for (j in seq_along(allow_params)) {
  for (i in seq_along(mean_shifts)) {
    N2 <- 0
    arl_track <- numeric(0)
    sd_track <- numeric(0)
    sd_arl <- sd_tol[i] + 1
    while (sd_arl > sd_tol[i]) {
      current_run_length <-
        f_out_control(
          lambda = allow_params[j],
          control_limit = cl_pred[j],
          OOC_dist_ops = list(mean = mean_shifts[i], sd =
                                1),
          tau=0,
          m=m,
          bootstrap_samples=3000
        )[["Time OOC"]]
      arl_track <- append(arl_track, current_run_length)
      sd_arl <- sd(arl_track) * 1.96 / sqrt(N2)
      arl_est <- mean(arl_track)
      if (is.na(sd_arl))
        sd_arl <- sd_tol[i] + 1
      N2 <- N2 + 1
      if (N2 < N2_min)
        sd_arl <- sd_tol[i] + 1
      if (N2 %% 40 == 0) {
        print(N2)
        print(arl_est)
        print(sd_arl)
      }
    }
    ooc_arl[i] <- arl_est
  }
  ARLs[[j]] <- ooc_arl
}


ooc_arl
saveresults("cvm_ewma_ooc_lambda")
# save.image("cvm_ewma_ooc_lambda.RData")

