# # Function archive
# power.res <- function(x, y, f, g = perm.test, boot = FALSE) {
#   # A nice function for power calculations Takes as input a two matrices. Each row a
#   # different draw of col's from the row dist.
#   dim.x <- dim(x)
#   fun <- f
#   pv <- NULL
#   if (boot) {
#     for (i in 1:dim.x[1]) {
#       a <- boot.test(x[i, ], y[i, ], f)
#       pv[i] <- a[[1]]
#     }
#     z <- pv
#     return(z)
#   }
#   for (i in 1:dim.x[1]) {
#     a <- g(x[i, ], y[i, ], f = fun)
#     pv[i] <- a[[1]]
#   }
#   z <- pv
#   return(z)
# }
# 
# perm.test <- function(x, y, ..., f, num.perm = 2001, diag = FALSE, exact = TRUE) {
#   # Runs a permutation test for a given function (Note: Input function should JUST
#   # output p-value) Takes as input f, a function that returns a statistic Computes the
#   # results of a permutation test for a test statistic Args: x: A vector of
#   # observations for a R.V. (must be numeric) y: Either (1) Another vector of
#   # observations (two sample) (2) A quantile function such as qnorm (one sample) f: a
#   # function which outputs a test statistic. MUST output only a numeric test statistic
#   # num.perm: Number of permutations. To avoid messiness, is best if not a multiple of
#   # 5 diag: if TRUE, also outputs the resulting values of the test statistics from
#   # each permutation exact: if TRUE, then for small samples <11 the function
#   # calculates every possible permutation Returns: A list containing: [1]: The p-value
#   # [2]: The critical value [3]: The test statistic [4](if diag=TRUE): The values of
#   # the test statistic in each permutation
#   require(gtools)
#   lenx <- length(x)
#   # This code snippet allows us to take in a function argument such as qgamma
#   if (is.character(y)) 
#     y <- get(y, mode = "function", envir = parent.frame())
#   if (is.function(y)) 
#     y <- y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...)  #Note: quantiles up for debate
#   leny <- length(y)
#   
#   # First step, combine into one dataset
#   z <- c(x, y)
#   lenz <- length(z)
#   # Calculate TS for the ACTUAL data
#   ts.obs <- f(x, y)
#   ts.random <- c(NULL)
#   ### 
#   if (lenz < 11 & exact == TRUE) {
#     all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, 
#                              set = FALSE)
#     all.permx <- all.perm[, 1:lenx]
#     all.permy <- all.perm[, (lenx + 1):lenz]
#     exact.perm <- dim(all.perm)[1]
#     for (i in 1:exact.perm) {
#       ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
#     }
#     p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
#     c.val <- quantile(ts.random, probs = 0.95)
#   } else {
#     for (i in 1:num.perm) {
#       z1 <- sample(z, size = lenz, replace = FALSE)
#       a <- z1[1:lenx]
#       b <- z1[(lenx + 1):lenz]
#       ts.random[i] <- f(a, b)
#     }
#     p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
#     c.val <- quantile(ts.random, probs = 0.95)
#   }
#   
#   # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
#   # statistic
#   if (diag == TRUE) {
#     return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
#   } else {
#     return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
#   }
# }
# 
# # 
# # wave.bec2 <- function(x,y, ..., interp = 4, doplot=F, wf="haar", reduce=2)
# # {
# #   library(wavethresh)
# #   library(waveslim)
# #   x <- sort(x)
# #   lenx <- length(x)
# #   # Two Sample Test
# #   if (is.numeric(y)) {
# #     leny <- length(y)
# #     y <- sort(y)
# #     
# #     num_quan <- min(2^floor(log(lenx/reduce,2)), 2^floor(log(leny/reduce,2)))
# #     prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
# #     
# #     # Because # of quantiles is < data points, need to interpolate
# #     qx <- quantile(x, probs = prob, type = interp)
# #     qy <- quantile(y, probs = prob, type = interp)
# #     q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
# #     q_inter_y <- approxfun(prob, qy, yleft = min(qy), yright = max(qy))
# #     quan_dif <- q_inter_y(prob) - q_inter_x(prob)
# #     print(quan_dif)
# #     n <- length(quan_dif)
# #     wave_tran <- wd(quan_dif, family="DaubExPhase", filter.number=1) #haar wavelet
# #     scoef <- sort(abs(wave_tran$D), decreasing=TRUE)
# #     sum_coef <- sum(abs(scoef))
# #     num_coef <- min(which(cumsum(abs(scoef))/sum_coef>.9))
# #     wave_th <- threshold(wave_tran, mannum= num_coef)
# #     STAT <- sum(wave_th$D^2)
# #     return(STAT)
# #   }
# #   
# #   # One Sample
# #   if (is.list(y)) 
# #     y <- names(y)
# #   if (is.function(y)) 
# #     funname <- as.character(substitute(y))
# #   if (is.character(y)) 
# #     funname <- y
# #   y <- get(funname, mode = "function", envir = parent.frame())
# #   if (!is.function(y)) 
# #     stop("'y' must be numeric or a function or a string naming a valid function")
# #   # Assuring diadic lengths
# #   #nx_p2 <- 2^floor(log(lenx,2))
# #   
# #   num_quan <- 2^floor(log(lenx,2))
# #   prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
# #   
# #   
# #   # Estimated Quantiles:
# #   qx <- quantile(x, probs = prob, type = interp)
# #   q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
# #   x_quans <- q_inter_x(prob)
# #   # True Quantiles
# #   true_quan <- y(prob, ...)
# #   
# #   quan_dif <- x_quans - true_quan
# #   print(quan_dif)
# #   
# #   n <- length(quan_dif)
# #   wave_tran <- wd(quan_dif, family="DaubExPhase", filter.number=1) #haar wavelet
# #   test_wave <- dwt(quan_dif,wf="haar")
# #   return(list(test_wave,wave_tran))
# #   
# #   print(quan_dif)
# #   scal_coef <- tail(wave_tran$C,1) # Scaling Coefficient
# #   print(scal_coef)
# #   # Take abs. value of wavelet coefficients, sort in decreasing order
# #   w_coef <- sort(abs(c(scal_coef,wave_tran$D)), decreasing=TRUE) 
# #   # Sum of them
# #   sum_coef <- sum(w_coef)
# #   # Get number of coefficients to keep, based on 90% thresholding
# #   num_coef <- min(which(cumsum(abs(w_coef))/sum_coef>.9))
# #   # Threshold according
# #   wave_th <- threshold(wave_tran, policy="mannum",value= num_coef)
# #   return(wave_th)
# #   # Calculate square of remaining wavelet coefficients
# #   th_scal_coef <- tail(wave_th$C,1)^2
# #   th_w_coef <- c(th_scal_coef, wave_th$D)
# #   print(th_w_coef)
# #   STAT <- sum(th_w_coef^2)
# #   return(STAT)
# # }

# wave.den.perm <- function (x=rnorm(10), y=rnorm(10), n=2^5, p=100, doplot=F, doplot1=F) 
# {
#   z <- c(x, y)
#   n.x <- length(x)
#   n.z <- length(z)
#   ts.vec <- matrix(,nrow=p, ncol=2)
#   
#   for(i in 1:p)
#   {
#     perm <- sample(z, n.z, replace=F)
#     #a <- wave.den(perm[1:n.x], perm[(n.x + 1):n.z])
#     #b <- ks.res.simp(perm[1:n.x], perm[(n.x + 1):n.z])
#     #n=n, 
#     #doplot=doplot) 
#     a <- wave.energy(perm[1:n.x], perm[(n.x + 1):n.z], n=n, 
#                      doplot=doplot) 
#     ts.vec[i, ] <- a
#   }
#   
#   colnames(ts.vec) <- c("oc", "ks")
#   rownames(ts.vec) <- rep("", p)
#   
#   d.oc <- density(ts.vec[, 1])
#   d.ks <- density(ts.vec[, 2])
#   
#   pp <- ceiling(0.95 * p)
#   oc.crit <- sort(ts.vec[, 1])[pp]
#   ks.crit <- sort(ts.vec[, 2])[pp]
#   # Plots
#   a <- wave.den(x, y, n=n, doplot=doplot)
#   
#   if(doplot1)
#   {
#     plot(d.oc, xlim=range(d.oc$x, d.ks$x), 
#          ylim=range(d.oc$y, d.ks$y),main="")
#     lines(d.ks, col=3)
#     abline(v=oc.crit, lty=3)
#     abline(v=ks.crit, col=3, lty=3)
#     abline(v=a[1])
#     abline(v=a[2], col=3)
#   }
#   
#   oc.dec <- ifelse(a[1] < oc.crit, 0, 1)
#   ks.dec <- ifelse(a[2] < ks.crit, 0, 1)
#   c(oc.dec, ks.dec)
# }
# 
# wave.den.power <- function (n=2^5, p=500, doplot=F, m=100, doplot1=F) 
# {
#   b <- numeric(0)
#   power.vec <- matrix(,nrow=m, ncol=2)
#   for(i in 1:m)
#   {
#     x <- rnorm(100)
#     #  	y <- rt(100, df=20)
#     y <- rnorm(100, sd=1)
#     a <- wave.den.perm(x=x, y=y, p=p)#, n=n, p=p, doplot=doplot, doplot1=doplot1)
#     power.vec[i, ] <- a 
#   }
#   c(oc=mean(power.vec[, 1]), ks=mean(power.vec[, 2]), frac=mean(power.vec[, 1]) / mean(power.vec[, 2]))
# }
# 
# boot.test <- function(x, y, f, num.perm = 1000, diag = FALSE, exact = TRUE) {
#   # Runs a BOOTSTRAP test for a given function (Note: Input function should JUST
#   # output p-value) Takes as input f, a function that returns a statistic
#   require(gtools)
#   lenx <- length(x)
#   leny <- length(y)
#   # First Step, combine into one dataset
#   z <- c(x, y)
#   lenz <- length(z)
#   # Calculate TS for the ACTUAL data
#   ts.obs <- f(x, y)
#   ts.random <- c(NULL)
#   ### 
#   if (lenz < 10 & exact == TRUE) {
#     all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, 
#                              set = FALSE)
#     all.permx <- all.perm[, 1:lenx]
#     all.permy <- all.perm[, (lenx + 1):lenz]
#     exact.perm <- dim(all.perm)[1]
#     for (i in 1:exact.perm) {
#       ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
#     }
#     p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
#     c.val <- quantile(ts.random, probs = 0.95)
#   } else {
#     for (i in 1:num.perm) {
#       z1 <- sample(z, size = lenz, replace = TRUE)
#       a <- z1[1:lenx]
#       b <- z1[(lenx + 1):lenz]
#       ts.random[i] <- f(a, b)
#     }
#     p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
#     c.val <- quantile(ts.random, probs = 0.95)
#   }
#   
#   # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
#   # statistic
#   if (diag == TRUE) 
#     return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random)) else {
#       return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
#     }
# }