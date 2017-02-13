# Testing the approximation of the CVM test for inverse gamma
library(statmod)
set.seed(5)
samp_size <- 50
x <- rnorm(samp_size)
STAT <- cvm.res(x, pnorm)
perm.test(x, pnorm, f=cvm.res)
mu1 <- 1/6
mu2 <- (4*samp_size-3)/ (180*samp_size)
mu2_corrected <- 5*samp_size /(24*samp_size-18)
pinvgamma(STAT, mu1, mu2)
pinvgamma(STAT, mu1, mu2_corrected)


# quick test
runs <- 100
test <- numeric(runs)
test2 <- numeric(runs)

for(i in 1:runs){
  samp_size <- 50
  x <- rnorm(samp_size)
  STAT <- cvm.res(x, pnorm)
  test2[i] <- perm.test(x, pnorm, f=cvm.res)[["p-value"]]
  mu1 <- 1/6
  mu2 <- (4*samp_size-3)/ (180*samp_size)
  mu2_corrected <- 5*samp_size /(24*samp_size-18)
  #pinvgamma(STAT, mu1, mu2)
  test[i] <- 1-pinvgauss(STAT, mu1, mu2_corrected)
}
par(mfrow=c(1,2))
plot(test, main="approximation",type="l")
points(test2, main="perm test", type="l", col="red")
hist(test)

### Anderson Darling Approximation
source("functions.R")
library(statmod)
set.seed(5)
samp_size <- 5
x <- rnorm(samp_size)

# quick test
runs <- 100
test <- numeric(runs)
test2 <- numeric(runs)

for(i in 1:runs){
  samp_size <- 5
  x <- rnorm(samp_size)
  STAT <- ad.res(x, pnorm)
  test2[i] <- perm.test(x, pnorm, f=ad.res)[["p-value"]]
  test[i] <- find_ad_pval_asym(STAT)
}
par(mfrow=c(1,2))
plot(test, main="approximation",type="l")
points(test2, main="perm test", type="l", col="red")
hist(test)