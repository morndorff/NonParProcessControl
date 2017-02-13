# Qiu's CUSUM chart
rm(list=ls())
set.seed(5)
source("functions.R")

start_time <- Sys.time()

# Setting Parameters
M <- 500
rdist <- rnorm
m <- 5
ARL <- 200
########

##### Parameters for UCL/LCL
# Most Important
allow_params = c(.25, .2, .15, .1, .05) 
ucl <- c(12, 12, 12, 12, 12)
lcl <- c(11.45, 12.5, 10.5, 10.5, 10.5)

num_allow_param <- length(allow_params)

# Tolerance Parameters
sd_tol <- 10
N2 <- 0 
N2_min <- 300 # At least 300 runs
sd_arl <- sd_tol + 2

arl_track <- numeric(0)
sd_track <- numeric(0)

cl_pred <- numeric(num_allow_param)
cl_int <- matrix(nrow=num_allow_param, ncol=3) #lower and upper confidence band
colnames(cl_int) <- c("fit", "lwr", "upr")

num_points <- 5 # For regression estimate
x <- matrix(nrow = num_allow_param, ncol = num_points)
y <- matrix(nrow = num_allow_param, ncol = num_points)
for(i in 1: num_allow_param){
  x[i , ] <- seq(lcl[i], ucl[i], length.out=num_points)
}
# Fill in x matrix

for(i in 1:num_allow_param){
  for(j in seq_along(x[i, ])){
    N2 <- 0 # At least 300 runs
    sd_arl <- sd_tol + 2
    arl_track <- NULL
    while(sd_arl > sd_tol){
      current_run_length <- qiu_ARL(control_limit=x[i, j], 
                                    kP=allow_params[i], m=5,
                                    ic_data=rnorm(M),
                                    num_bps=5)[["Time OOC"]]
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
      y[i, j] <- mean(arl_track)
    }
  }
  xi <- x[i,] 
  yi <- y[i,]
  
  m1 <- lm(xi ~ yi + I(yi^2))
  cl_pred[i] <- predict(m1, newdata=data.frame(yi=ARL))
  paste("For allowance parameter:", allow_params, "Calculated ARL200 is", cl_pred[i])
  
  cl_int[i, ] <- predict(m1, newdata=data.frame(yi=ARL), interval="confidence")
  print(paste("X-Values are", x[i, ]))
  print(paste("Y-Values are", y[i, ]))
}




############## could wrap another loop around this for different values of the allowance parameter
liAllowParams_And_ControlLimits <- list()




mean_shifts <- c(1, .9, .8, .7 , .6, .5, .4, .3, .2)
sd_tol <- 1:10 # temporary
ooc_arl <- numeric(0)
ARLs <- list(num_allow_param)
for(j in seq_along(allow_params)){
  
  for(i in seq_along(mean_shifts)){
    N2 <- 0
    arl_track <- numeric(0)
    sd_track <- numeric(0)
    sd_arl <- sd_tol[i] + 1
    while(sd_arl > sd_tol[i]){
      current_run_length <- qiu_Phase_II(kP=allow_params[j], 
                                         control_limit=cl_pred[j], 
                                         OOC_dist_ops=list(mean=mean_shifts[i], 
                                                           sd=1))[["Time OOC"]]
      arl_track <- append(arl_track, current_run_length)
      sd_arl <- sd(arl_track)* 1.96 / sqrt(N2)
      arl_est <- mean(arl_track)
      if(is.na(sd_arl)) sd_arl <- sd_tol[i] + 1
      N2 <- N2 + 1
      if(N2 < N2_min) sd_arl <- sd_tol[i] + 1
      if(N2 %% 40==0){
        print(N2)
        print(arl_est)
        print(sd_arl)
      }
    }
    ooc_arl[i] <- arl_est
  }
  ARLs[[j]] <- ooc_arl
}

names(ooc_arl) <- mean_shifts
plot(mean_shifts, log(ooc_arl), main="Allowance Parameter = 2")
plot(mean_shifts, ooc_arl, main="Allowance Parameter = 2")

save.image("arl200_new.RData")


















# 
# 
# # Establish initial range
# for(i in seq_along(allowance_parameter)){
#   # 3 Stabs with regression
#   for(j in 1:num_iter){
#     res <-
#       Find_ARL_MC_Reg(
#         f = qiu_CUSUM_generic,
#         tstat=SP_Two_Quick,
#         allowance_param = allowance_parameter[i],
#         m = 5,
#         Cores = Cores,
#         ucl = ucl_vec[i],
#         lcl = lcl_vec[i],
#         arl = ARL,
#         num_points = num_points[j],
#         sd_tol = ARL / tolerance[j],
#         bootstrap_samples=2500
#       )
#     print(res)
#     point_pred <- res[["Prediction"]]
#     interval <- res[["Confidence Interval"]]
#     lcl_vec[i] <- interval[2]
#     ucl_vec[i] <- interval[3]
#   }
#   calculated_CL[i] <- point_pred
# }