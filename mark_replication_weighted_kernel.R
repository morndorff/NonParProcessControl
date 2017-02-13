# Kernel based method -- gives pretty much the same results as the weighted method, which isn't that suprising.

rm(list=ls())
set.seed(5)
source("functions.R")

library(doMC)
library(parallel)

# Params 
M <- 500 # In control sample size
tstat <- "cvm"

dist_tested <- "normal"
rdist <- standard_normal
rdist_character = "standard_normal"
ic_params <- list(mean=0, sd=1)

ARL <- 200
m <- 5 # sample size at each iteration

# Parallel Processing Parameters
Cores <- 3
registerDoMC(Cores)

# Find Control Control Limits for the Process
# Most Important
f_in_control <- Find_CL_CVM_Weight_Kernel
f_out_control <- Find_CL_CVM_Weight_Kernel_OOC

# TEMPORARY
method <- "weighted_kernel"
lambdas=c(.01)
# lambdas = c(.01, .02, .03) # double window change
num_lambdas <- length(lambdas) 

ucl = c(.04)
lcl = c(.006)
# ucl = c(.04, .04, .04)
# lcl = c(.01, .01, .01)

sd_tol_init = 10
data.results <- NULL

# Initial Control Limit Finding -------------------------------------------
#

find_cl <- function(lambdas=lambdas, ucl, lcl, rdist, rdist_character, ic_params, sd_tol_init, method, ...){
  
  # Initial Estimates for Lambdas
  #ucl <- c(.02, .02, .02, .02, .02)
  names(ucl) <- as.character(lambdas)
  #lcl <- c(.005, .005, 0.004, 0.004, 0.004)
  names(lcl) <- as.character(lambdas)
  
  CL_200_Final <- numeric(num_lambdas)
  names(CL_200_Final) <- as.character(lambdas)
  # Its opposite now!
  for(i in seq_along(lambdas)){
    control_limit_sequence <- seq(ucl[as.character(lambdas[i])], 
                                  lcl[as.character(lambdas[i])], 
                                  length.out=50)
    
    # TODO: add code to do different distribution!
    cl.estimate <- CL_Finding_With_UB(allow_param= lambdas[i], 
                                      cl_seq=control_limit_sequence, 
                                      f_in_control=f_in_control, 
                                      sd_tol_init=sd_tol_init, 
                                      Cores=3,
                                      max_iter=5000,
                                      pval_approach=TRUE,
                                      method=method,
                                      ICdist=rdist_character,
                                      IC_dist_ops=ic_params)
    cl.estimate.avg <- rowMeans(cl.estimate)
    CLs <- as.numeric(names(cl.estimate.avg))
    f <- approxfun(cl.estimate.avg, CLs)
    plot(cl.estimate.avg, CLs)
    plot(CLs, cl.estimate.avg)
    CL_200_Final[i] <- f(ARL)
  }
  CL_200_Final
}

# Out of control function -------------------------------------------------


find_ooc_limits <-  function(shifts,
                             sd_tol = 3,
                             rdist. = rdist, #error catch
                             rIC = rdist_character,
                             IC_dist_ops= ic_params,
                             rOOC = rdist_character,
                             type="mean",
                             N2min=3000,
                             N2_max=4000,
                             lambdas,
                             climit,
                             method=method,
                             output=NULL)  {
  # Function documentation
  # shifts -- the amount of mean shift -- for now
  # sd_tol -- how good should the estimate should be -- right now unnecessary, see n2min/n2max
  # Distribution options -- how should things move
  # CL_200_Final -- a vector of control limit
  # lambdas -- a vector of parameters
  # output: put results out as a data frame? if so, output="dataframe"
  sd_tol <- rep(sd_tol, length(shifts))
  # sd_tol <- sd_tol * 2 #temporary
  num_lambdas <- length(lambdas)
  num_alternatives <- length(shifts)
  
  ooc_arl <- numeric(0)
  ooc_sd <- numeric(0)
  ARLs <- list(num_lambdas)
  ARL_sd <- list(num_lambdas)
  
  N2_min = 3000 # replace later
  Num_Iter <- 6 * Cores
  
  for (j in seq_along(lambdas)) {
    for (i in seq_along(shifts)) {
      
      if(type=="mean"){
        OOC_dist_ops <- c(append(IC_dist_ops, list(center_mean=shifts[i])))
      }else if(type=="sd"){
        OOC_dist_ops <- c(append(IC_dist_ops, list(center_sd=shifts[i])))
      }
      
      N2 <- 0
      arl_track <- numeric(0)
      sd_track <- numeric(0)
      sd_arl <- sd_tol[i] + 1
      while (sd_arl > sd_tol[i]) {
        
        par_obj <- foreach(n=1:Num_Iter) %dopar%{
          current_run_length <-
            f_out_control(
              ic_data= do.call(rdist., c(list(M), IC_dist_ops)),
              lambda = lambdas[j],
              control_limit = climit[j],
              ICdist = rIC,
              OOCdist = rOOC,
              IC_dist_ops = IC_dist_ops,
              OOC_dist_ops = OOC_dist_ops,
              tau=0,
              m=m,
              method=method
            )[["Time OOC"]] 
        }
        par_obj <- unlist(par_obj)
        
        arl_track <- append(arl_track, par_obj)
        N2 <- N2 + Num_Iter #MC
        sd_arl <- sd(arl_track) * 1.96 / sqrt(N2)
        arl_est <- mean(arl_track)
        if (is.na(sd_arl))sd_arl <- sd_tol[i] + 1 
        if (N2 < N2_min)  sd_arl <- sd_tol[i] + 1
        if (N2 > N2_max)  sd_arl <- sd_tol[i] - 1
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
  
  if(output=="dataframe"){
    # attributes to keep track of: lambda, ooc_type, shifts, sd, ic_distribution, parameters
    out_data = NULL
    for(i in 1:num_lambdas){
      ooc_time = results_matrix[i,]
      temp_data <- as.data.frame(ooc_time)
      temp_data$sd = results_matrix_sd[i, ]
      temp_data$lambda = lambdas[i]
      temp_data$ooc_type = ooc_type
      temp_data$shifts = shifts
      temp_data$ic_dist = rdist_character
      temp_data$param = paste(as.character(ic_params), collapse="", sep="")
      temp_data$batchsize = m
      temp_data$icsize = M
      temp_data$ARL= ARL
      temp_data$method= method
      if(is.null(out_data)){
        out_data = temp_data
      }else{
        out_data = rbind(out_data, temp_data)
      }
      temp_data = NULL
    }
    return(out_data)
  }
  return(list("Results Matrix"=results_matrix, "SD Matrix"=results_matrix_sd))
  
}

# Graph Function ----------------------------------------------------------



generate_graph <- function(results_matrix= results[[1]], lambdas=lambdas, shift=shifts, title_text){
  library(tidyverse)
  library(ggplot2)
  plot_data <- as.data.frame(cbind(lambdas, results_matrix))
  plot_data <- gather(plot_data, shift, ARL, -lambdas)
  plot_data$lambdas <- as.factor(plot_data$lambdas)
  plot_data$shift <- as.numeric(plot_data$shift)
  ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
    geom_line() +
    geom_point() +
    coord_trans(y = "log") +
    
    if(!is.null(title_text)) ggtitle(title_text)
}




# Find the control limit --------------------------------------------------

control_limit <-
  find_cl(
    lambdas = lambdas,
    ucl = ucl,
    lcl = lcl,
    rdist = rdist,
    rdist_character = rdist_character,
    ic_params = ic_params,
    method = method,
    sd_tol_init=sd_tol_init)


# OOC Mean Shifts ----------------------------------------------------

ooc_type <- "mean"
shifts <- seq(-1, 1, by=.2)
shifts <- c(-1, -.8, -.6, -.4, -.2)

results <- find_ooc_limits(
           rdist. = rdist,
           shifts = shifts,
           sd_tol=3,
           rOOC = rdist_character,
           rIC = rdist_character,
           IC_dist_ops = ic_params,
           type=ooc_type,
           N2min=3000,
           N2_max=4000,
           lambdas= lambdas,
           climit=control_limit,
           method=method,
           output="dataframe")

# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)



# Generate Graph ----------------------------------------------------------
# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )


# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)



# SD SHIFTS ---------------------------------------------------------------
ooc_type <- "sd"
shifts <- seq(.2, 2, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)

# Generate Graph ----------------------------------------------------------


# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------



# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)









# Find Control Limit ------------------------------------------------------

#
# chi_sq(1) distribution
#

dist_tested <- "chi_square"
rdist <- standard_chi
rdist_character = "standard_chi"
ic_params <- list(df=1)



control_limit <-
  find_cl(
    lambdas = lambdas,
    ucl = ucl,
    lcl = lcl,
    rdist = rdist,
    rdist_character = rdist_character,
    ic_params = ic_params,
    method = method,
    sd_tol_init=sd_tol_init)



# Mean Shift --------------------------------------------------------------


ooc_type <- "mean"
shifts <- seq(-1, 1, by=.2)


results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)

# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)



# SD SHIFTS ---------------------------------------------------------------
ooc_type <- "sd"
shifts <- seq(.2, 2, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)

# Generate Graph ----------------------------------------------------------


# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)







# Find Control Limit ------------------------------------------------------

#
# t(4) distribution
#

dist_tested <- "t"
rdist <- standard_t
rdist_character = "standard_t"
ic_params <- list(df=4)




control_limit <-
  find_cl(
    lambdas = lambdas,
    ucl = ucl,
    lcl = lcl,
    rdist = rdist,
    rdist_character = rdist_character,
    ic_params = ic_params,
    method = method,
    sd_tol_init=sd_tol_init)



# Mean Shift --------------------------------------------------------------


ooc_type <- "mean"
shifts <- seq(-1, 1, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)


# Generate Graph ----------------------------------------------------------

# 
# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)


# SD SHIFTS ---------------------------------------------------------------
ooc_type <- "sd"
shifts <- seq(.2, 2, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)

# Generate Graph ----------------------------------------------------------


# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)












# Find Control Limit ------------------------------------------------------
#
# chisq(4) distribution
#

dist_tested <- "chisq"
rdist <- standard_chi
rdist_character = "standard_chi"
ic_params <- list(df=4)


control_limit <-
  find_cl(
    lambdas = lambdas,
    ucl = ucl,
    lcl = lcl,
    rdist = rdist,
    rdist_character = rdist_character,
    ic_params = ic_params,
    method = method,
    sd_tol_init=sd_tol_init)



# Mean Shift --------------------------------------------------------------


ooc_type <- "mean"
shifts <- seq(-1, 1, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)


# Generate Graph ----------------------------------------------------------


# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)


# SD SHIFTS ---------------------------------------------------------------
ooc_type <- "sd"
shifts <- seq(.2, 2, by=.2)

results <- find_ooc_limits(
  rdist. = rdist,
  shifts = shifts,
  sd_tol=3,
  rOOC = rdist_character,
  rIC = rdist_character,
  IC_dist_ops = ic_params,
  type=ooc_type,
  N2min=3000,
  N2_max=4000,
  lambdas= lambdas,
  climit=control_limit,
  method=method,
  output="dataframe")
# Generate data frame
# Merge results with others
data.results = rbind(data.results, results)

# Generate Graph ----------------------------------------------------------

# 
# generate_graph(
#   results_matrix = results[[1]],
#   lambdas = lambdas,
#   shift = shifts,
#   title_text=paste("ARL Values for",ooc_type,"Shift,",method," Method;", "In Control Distribution:", dist_tested)
# )

# Save Results ------------------------------------------------------------

# results_name = paste(ooc_type,"_",tstat,"_",method, "_arl_",200,"_icdist_",dist_tested,"param",paste(as.character(ic_params), collapse="", sep=""), sep="")
# 
# README_OOC <- paste(ooc_type, "SHIFT RESULTS FOR CVM PVAL,",method, "METHOD, ARL 200, distribution tested", dist_tested, "parameter values:", paste(as.character(ic_params), collapse="", sep=""))
# save_name <-nameresults(results_name)
# save(control_limit, lambdas, README_OOC, ARL, m, M, sd_tol_init, 
#      shifts, results, method,dist_tested,  file=save_name)



# Save Final Results ------------------------------------------------------


results_name = paste(tstat,"_",method, "_arl_",ARL, "batch", m,sep="")

save_name <-nameresults(results_name)
save(data.results, tstat, method, ARL, m, M, file=save_name)

