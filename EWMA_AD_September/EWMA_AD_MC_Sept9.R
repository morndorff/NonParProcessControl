# This script identifies the control limits for a Cramer Von Mises based 
# EWMA control chart. It does so for various values of lambda.

# It then additionally identifies the out of control average run lengths

rm(list=ls())
setwd("/home/mark/Dropbox/Research/Method_Sims/")
source("functions.R")
set.seed(5000)
library(doMC)
library(foreach)

README_CL <- "CONTROL LIMITS FOR A-D EWMA CONTROL CHART"


# Params 
M <- 500 # In control sample size
rdist <- rnorm
rdist_character <- "rnorm"
ic_params <- list(mean=0, sd=1)

# Adjusted Distribution
adist <- rt
adist_character <- "rt"
a_params = list(df=5)

ARL <- 200
m <- 5 # sample size at each iteration

# Parallel Processing Parameters
Cores <- 3
registerDoMC(Cores)

# Find Control Control Limits for the Process
##### Parameters for UCL/LCL
# Most Important
f_in_control <- EWMA_Find_CL_AD
f_out_control <- EWMA_AD_OOC

# Initial Estimates for Lambdas
lambdas = c(.05, .1, .2, .3)
num_lambdas <- length(lambdas) 
# from df_adjust=5
# init.lcl=c(.3, .45, .75, 1.02)
# init.ucl=c(.45, .7, 1.15, 1.6)

init.lcl=c(0, 0, 0, 0)
init.ucl=c(.36, .61, 1.07, 1.5)



# Accuracy of estimate of UCLs
sd_tol_init = 7

# Writing new function for multiple allowance parameters


newccfun <- function(num_candidates=100, 
                     lcl, 
                     ucl, 
                     bootstrap_samples=6000, 
                     ic_sample_size=500,
                     allowance_parameter=.05,
                     num_replications=10,
                     target_ARL=200,
                     runs=1000,
                     df_opt=TRUE,
                     df_adjust=5){
  ucl_mat <-matrix(, nrow=length(allowance_parameter), ncol=3)
  cl_mat <- matrix(nrow=2, ncol=length(lcl))
  cl_mat[1, ] = lcl
  cl_mat[2, ] = ucl
  rownames(cl_mat) <- c("lcl", "ucl")
  for(num in seq_along(allowance_parameter)){
    estimate_ucl_200 <- numeric(num_replications)
    estimate_ucl_100 <- numeric(num_replications)
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
        in_control_tstat <- In_Control_Dist_AD_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
                                                    ICdist="rt", ic_params=list(df=df_adjust*allowance_parameter[num] / .07)) ## HARD CODE CITY, CAPTAIN 
      }else{
      in_control_tstat <- In_Control_Dist_AD_Exp(ic_data=ic_data, m=5, bootstrap_samples=bootstrap_samples,
                                                  ICdist="rt", ic_params=list(df=df_adjust))
      }
      # get mean of the distribution
      mean_ic_tstat <- mean(in_control_tstat)
      print(paste("mean tstat",mean_ic_tstat))
      # hist(in_control_tstat)
      #plot(density(in_control_tstat))
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
      estimate_ucl_100[reps] <- f_in(100)
      #ucl_mat[num, 2] <- estimate_ucl_180
      estimate_ucl_250[reps] <- f_in(250)
      #ucl_mat[num, 3] <- estimate_ucl_250
    }
    ucl_mat[num, 1] <- mean(estimate_ucl_200)
    ucl_mat[num, 2] <- mean(estimate_ucl_100)
    ucl_mat[num, 3] <- mean(estimate_ucl_250)
    
  }
  colnames(ucl_mat) <- c("200", "100", "250")
  rownames(ucl_mat) <- as.character(allowance_parameter)
  ucl_mat
}

num_replications <- 10
control_limit_mat <- newccfun(ic_sample_size=M, allowance_parameter = lambdas, lcl=init.lcl, ucl=init.ucl,
                              num_replications = num_replications, df_adjust=20)
control_limit_mat
is.na(control_limit_mat)

library(beepr)
beep()



lcl <- control_limit_mat[, "100"]
ucl <- control_limit_mat[, "250"]

# Check to see if the UCL for 250 is accurate
#If it is, use a grid approach to search for ARL 200

CL_200_Final <- numeric(num_lambdas)
names(CL_200_Final) <- as.character(lambdas)
for(i in seq_along(lambdas)){
  control_limit_sequence <- seq(ucl[as.character(lambdas[i])], 
                                lcl[as.character(lambdas[i])], 
                                length.out=50)
  cl.estimate <- CL_Finding_With_UB(allow_param= lambdas[i], 
                     cl_seq=control_limit_sequence, 
                     f_in_control=f_in_control, 
                     sd_tol_init=sd_tol_init, 
                     Cores=3,
                     max_iter=5000)
  cl.estimate.avg <- rowMeans(cl.estimate)
  CLs <- as.numeric(names(cl.estimate.avg))
  f <- approxfun(cl.estimate.avg, CLs)
  plot(cl.estimate.avg, CLs)
  CL_200_Final[i] <- f(ARL)
}
beep()
save_name <-nameresults("ad_ewma_200")
save(CL_200_Final, lambdas, README_CL, ARL, m, M, sd_tol_init, file=save_name)

#######################
# BEGIN OOC SIMULATIONS
#######################


mean_shifts <- c(1, .9, .8, .7 , .6, .5, .4, .3, .2)
sd_tol <- c(.25, .3, .4, .45, .5, 2, 4, 6, 8)  # temporary
# sd_tol <- sd_tol * 2 #temporary

num_alternatives <- length(mean_shifts)

ooc_arl <- numeric(0)
ooc_sd <- numeric(0)
ARLs <- list(num_lambdas)
ARL_sd <- list(num_lambdas)

N2_min = 3000 # replace later
Num_Iter <- 6 * Cores

for (j in seq_along(lambdas)) {
  for (i in seq_along(mean_shifts)) {
    N2 <- 0
    arl_track <- numeric(0)
    sd_track <- numeric(0)
    sd_arl <- sd_tol[i] + 1
    while (sd_arl > sd_tol[i]) {
      par_obj <- foreach(n=1:Num_Iter) %dopar%{
        current_run_length <-
          f_out_control(
            lambda = lambdas[j],
            control_limit = CL_200_Final[j],
            OOC_dist_ops = list(mean = mean_shifts[i], sd = 1),
            tau=0,
            m=m,
            bootstrap_samples=6000
          )[["Time OOC"]]
      }
      par_obj <- unlist(par_obj)

      arl_track <- append(arl_track, par_obj)
      N2 <- N2 + Num_Iter #MC
      sd_arl <- sd(arl_track) * 1.96 / sqrt(N2)
      arl_est <- mean(arl_track)
      if (is.na(sd_arl))
        sd_arl <- sd_tol[i] + 1 

      if (N2 < N2_min)
        sd_arl <- sd_tol[i] + 1
      if (N2 %% 40 == 0) {
        print(N2)
        print(arl_est)
        print(sd_arl)
      }
    }
    ooc_arl[i] <- arl_est
    ooc_sd[i] <- sd_arl
    print(paste("Shift", mean_shifts[i], "done for lambda", lambdas[j], "which is", j, "of", length(lambdas)))
  }
  ARLs[[j]] <- ooc_arl
  ARL_sd[[j]] <- ooc_sd
}

results_matrix <- matrix(unlist(ARLs), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix) <- as.character(lambdas)
colnames(results_matrix) <- as.character(mean_shifts)

results_matrix_sd <- matrix(unlist(ARL_sd), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix) <- as.character(lambdas)
colnames(results_matrix) <- as.character(mean_shifts)

library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")


saveresults("ad_ewma_ooc_lambda")

README_OOC <- "MEAN SHIFT RESULTS FOR AD EWMA, ARL 200"
save_name <-nameresults("ad_ewma_200_mean_shift")
save(CL_200_Final, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
     mean_shifts, results_matrix, results_matrix_sd, N2_min, file=save_name)

#### SD Shifts ####
shifts <- c(-.5, -.25, -.1, .2, .35, .5, .75, 1)
sd_tol <- rep(8, length(shifts))  # temporary
N2_min = 3000 # replace later
N2_max = 2 # trouble shooting
# sd_tol <- sd_tol * 2 #temporary

num_alternatives <- length(shifts)

ooc_arl <- numeric(0)
ooc_sd <- numeric(0)
ARLs <- list(num_lambdas)
ARL_sd <- list(num_lambdas)


Num_Iter <- 6 * Cores

for (j in seq_along(lambdas)) {
  for (i in seq_along(shifts)) {
    N2 <- 0
    arl_track <- numeric(0)
    sd_track <- numeric(0)
    sd_arl <- sd_tol[i] + 1
    while (sd_arl > sd_tol[i]) {
      par_obj <- foreach(n=1:Num_Iter) %dopar%{
        current_run_length <-
          f_out_control(
            lambda = lambdas[j],
            control_limit = CL_200_Final[j],
            OOC_dist_ops = list(mean = 0, sd = 1 + shifts[i]),
            tau=0,
            m=m,
            bootstrap_samples=6000
          )[["Time OOC"]]
      }
      par_obj <- unlist(par_obj)
      
      arl_track <- append(arl_track, par_obj)
      N2 <- N2 + Num_Iter #MC
      sd_arl <- sd(arl_track) * 1.96 / sqrt(N2)
      arl_est <- mean(arl_track)
      if (is.na(sd_arl))
        sd_arl <- sd_tol[i] + 1 
      
      if (N2 < N2_min)  sd_arl <- sd_tol[i] + 1
      if(N2 > N2_max) sd_arl <- sd_tol[i] -1
      if (N2 %% 40 == 0) {
        print(N2)
        print(arl_est)
        print(sd_arl)
      }
    }
    ooc_arl[i] <- arl_est
    ooc_sd[i] <- sd_arl
    print(paste("Shift", shifts[i], "done for lambda", lambdas[j], "which is", j, "of", length(lambdas)))
  }
  ARLs[[j]] <- ooc_arl
  ARL_sd[[j]] <- ooc_sd
}

results_matrix <- matrix(unlist(ARLs), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix) <- as.character(lambdas)
colnames(results_matrix) <- as.character(shifts)

results_matrix_sd <- matrix(unlist(ARL_sd), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix_sd) <- as.character(lambdas)
colnames(results_matrix_sd) <- as.character(shifts)

cols <- rainbow(4)
plot(mean_shifts, results_matrix[1,], type="l", col=cols[1])
for(i in 2:dim(results_matrix)[1]){
  lines(mean_shifts, results_matrix[i,], type="l", col=cols[i])
}



library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line(size=2, aes(linetype=lambdas)) +
  coord_trans(y = "log")

README_OOC <- "SD SHIFT RESULTS FOR AD EWMA, ARL 200"
save_name <-nameresults("ad_ewma_200_sd_shift")
save(CL_200_Final, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
     shifts, results_matrix, results_matrix_sd, N2_min, file=save_name)
