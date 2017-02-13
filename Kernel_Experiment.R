# Experimenting with kernels for weighted data

set.seed(12)
x <- rlnorm(50, meanlog=1, sdlog=1)
## test if x follows a Gamma distribution with shape .6 and rate .1
dgeometric.test(x, dgamma, par=list(shape=0.6, rate=0.1), lower=0, upper=Inf, n.sim=100)
f0 <- function(x) ifelse(x>=0 & x<=1, 2-2*x, 0)
## test if risk76.1929 follows the distribution characterized by f0
dgeometric.test(risk76.1929, f0, lower=0, upper=1, n.sim=31)



fan.test(rbeta(100,1, 2), dunif, lower=0, upper=1)

f0 <- function(x) ifelse(x>=0 & x<=1, 2-2*x, 0)
## testing if risk76.1929 follows the distribution characterized by f0
fan.test(risk76.1929, f0, lower=0, upper=1, kernel="epanech")

set.seed(5)
source("functions.R")

x = rnorm(20)
lambda = .02
w = (1-lambda)^(1:20)
w= w/sum(w)
cvm.res.weight.internal.kernel(x, pnorm, w=w)
library(spatstat)
cvm.res.weight.internal(x, pnorm, w=w)

set.seed(5)
test = numeric(20)
for(i in 1:20){
ic_data <- rnorm(500)
test[i] = Fast_Bootstrap_CVM(ic_data= ic_data, m=5, bootstrap_samples=5000)
}
test
mean(test)
sd(test)


set.seed(5)
test = numeric(20)
for(i in 1:20){
  ic_data <- rnorm(500)
  test[i] = Fast_Bootstrap_CVM(ic_data= ic_data, m=5, bootstrap_samples=10000)
}
test
mean(test)
sd(test)


# expected change in standard deviation if mean of all entries is the same

set.seed(5)
test = numeric(20)
for(i in 1:20){
  ic_data <- rnorm(500)
  test[i] = Fast_Bootstrap_CVM(ic_data= ic_data, m=10, bootstrap_samples=5000)
}
test
mean(test)
sd(test)


set.seed(5)
test = numeric(20)
for(i in 1:20){
  ic_data <- rnorm(500)
  test[i] = Fast_Bootstrap_CVM(ic_data= ic_data, m=10, bootstrap_samples=10000)
}
test
mean(test)
sd(test)