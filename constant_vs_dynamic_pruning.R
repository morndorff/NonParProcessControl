# Constant vs Dynamic Pruning

# Tracking Mean Shift
rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_ewma_200_mean_shift_2016-10-17 143637.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")

rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_constant_ewma_200_mean_shift_2016-10-20 085549.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")

# Tracking SD Shift
rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_ewma_200_sd_shift_2016-10-17 144831.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
plot_data$shift <- as.numeric(plot_data$shift)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")

rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_constant_ewma_200_sd_shift_2016-10-20 091148.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data$lambdas <- as.factor(plot_data$lambdas)
plot_data$shift <- as.numeric(plot_data$shift)
ggplot(data=plot_data, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")


# Look at .03 keep

rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_constant_ewma_200_mean_shift_2016-10-20 085549.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data2 <- subset(plot_data, lambdas==.03)



load("~/Dropbox/Research/Method_Sims/cvm_pval_ewma_200_mean_shift_2016-10-17 143637.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)


plot_data_compare <- rbind(plot_data, plot_data2)

plot_data_compare$lambdas <- as.factor(plot_data_compare$lambdas)

ggplot(data=plot_data_compare, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log")+
  ggtitle("OOC ARL's for Mean Shift")




########## SD Compare
rm(list=ls())
load("~/Dropbox/Research/Method_Sims/cvm_pval_constant_ewma_200_sd_shift_2016-10-20 091148.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)
plot_data2 <- subset(plot_data, lambdas==.03)



load("~/Dropbox/Research/Method_Sims/cvm_pval_ewma_200_sd_shift_2016-10-17 144831.RData")
library(tidyverse)
library(ggplot2)
plot_data <- as.data.frame(cbind(lambdas, results_matrix))
plot_data <- gather(plot_data, shift, ARL, -lambdas)


plot_data_compare <- rbind(plot_data, plot_data2)

plot_data_compare$lambdas <- as.factor(plot_data_compare$lambdas)
plot_data_compare$shift <- as.numeric(plot_data_compare$shift)
ggplot(data=plot_data_compare, aes(shift, ARL, group=lambdas, color=lambdas)) +
  geom_line() +
  coord_trans(y = "log") +
  ggtitle("OOC ARL's for SD Shift")



