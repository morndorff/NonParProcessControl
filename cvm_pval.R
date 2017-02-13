# Test cvm p-value approximation process control
set.seed(5)
source("functions.R")

####### Finding an upper bound with this function
library(doMC)
library(parallel)
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
f_in_control <- EWMA_Find_CL_CVM_Pval
f_out_control <- EWMA_Find_CL_CVM_Pval_OOC

# Initial Estimates for Lambdas
lambdas = c(1,2,3,4,5)
num_lambdas <- length(lambdas) 
method = "partha" # Use partha prunign method

# Accuracy of estimate of UCLs
sd_tol_init = 3
ucl <- c(.02, .02, .02, .02, .02)
names(ucl) <- as.character(lambdas)
lcl <- c(.005, .005, 0.004, 0.004, 0.004)
names(lcl) <- as.character(lambdas)


CL_200_Final <- numeric(num_lambdas)
names(CL_200_Final) <- as.character(lambdas)
# Its opposite now!
for(i in seq_along(lambdas)){
  control_limit_sequence <- seq(ucl[as.character(lambdas[i])], 
                                lcl[as.character(lambdas[i])], 
                                length.out=50)
  cl.estimate <- CL_Finding_With_UB(allow_param= lambdas[i], 
                                    cl_seq=control_limit_sequence, 
                                    f_in_control=f_in_control, 
                                    sd_tol_init=sd_tol_init, 
                                    Cores=3,
                                    max_iter=5000,
                                    pval_approach=TRUE)
  cl.estimate.avg <- rowMeans(cl.estimate)
  CLs <- as.numeric(names(cl.estimate.avg))
  f <- approxfun(cl.estimate.avg, CLs)
  plot(cl.estimate.avg, CLs)
  plot(CLs, cl.estimate.avg)
  CL_200_Final[i] <- f(ARL)
}




### OOC Simulations

shifts <- c(1, .9, .8, .7 , .6, .5, .4, .3, .2)
sd_tol <- c(.25, .3, .4, .45, .5, 2, 4, 6, 8)  # temporary
# sd_tol <- sd_tol * 2 #temporary

num_alternatives <- length(shifts)

ooc_arl <- numeric(0)
ooc_sd <- numeric(0)
ARLs <- list(num_lambdas)
ARL_sd <- list(num_lambdas)

N2_min = 3000 # replace later
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
            OOC_dist_ops = list(mean = shifts[i], sd = 1),
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
    print(paste("Shift", shifts[i], "done for lambda", lambdas[j], "which is", j, "of", length(lambdas)))
  }
  ARLs[[j]] <- ooc_arl
  ARL_sd[[j]] <- ooc_sd
}

results_matrix <- matrix(unlist(ARLs), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix) <- as.character(lambdas)
colnames(results_matrix) <- as.character(shifts)

results_matrix_sd <- matrix(unlist(ARL_sd), nrow=num_lambdas, ncol= num_alternatives, byrow=TRUE)
rownames(results_matrix) <- as.character(lambdas)
colnames(results_matrix) <- as.character(shifts)

library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")

README_OOC <- "Mean SHIFT RESULTS FOR CVM PVAL, PARTHA METHOD, ARL 200"
save_name <-nameresults("cvm_pval_ewma_200_mean_shift")
save(CL_200_Final, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
     shifts, results_matrix, results_matrix_sd, N2_min, file=save_name)




########## SD Shifts

#### SD Shifts ####
shifts <- seq(-.8, 1, by=.2)
shifts <- shifts[-which(shifts==0)]
sd_tol <- rep(8, length(shifts))  # temporary
N2_min = 3000 # replace later
N2_max = 4000
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
plot(shifts, results_matrix[1,], type="l", col=cols[1])
for(i in 2:dim(results_matrix)[1]){
  lines(shifts, results_matrix[i,], type="l", col=cols[i])
}



library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line(size=2, aes(linetype=lambdas)) +
  coord_trans(y = "log")

README_OOC <- "SD SHIFT RESULTS FOR CVM PVAL, PARTHA METHOD, ARL 200"
save_name <-nameresults("cvm_pval_ewma_200_sd_shift")
save(CL_200_Final, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
     shifts, results_matrix, results_matrix_sd, N2_min, file=save_name)

