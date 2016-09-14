# Useful Functions
source("functions_out.R")
source("functions_plot.R")
source("functions_power.R")
source("functions_string.R")
source("functions_tstats.R")
source("functions_wavelets.R")
source("functions_SPS.R")
source("functions_MC.R")
source("functions_fourier.R")
source("functions_cv.R")
source("functions_sp.R")
source("functions_uniform.R")
source("functions_qiu.R")
source("functions_package.R")
source("functions_ewma.R")

perm.test <- function(x, y, distops = NULL, f, fops = NULL, num.perm = 1501, diag = FALSE, 
                      exact = FALSE, out=FALSE, do.plot=FALSE, ...) {
  #Args: 
  # x: numeric vector
  # y: numeric vector or quantile distribution (qgamma, qnorm)
  # distops: list describing quantile parameters (e.g. for qunif, list(min=0,max=2)
  # f: function outputting a test statistic
  # fops: if the test statistic has options, put them here. (e.g. for myts.out, size=.2)
  # num.perm: number of permutations to assess p-values
  # 
  # Output:
  # list containing observed test statistic, 
  #if (is.null(distops)==FALSE){
  #  if(is.list(distops)==FALSE) stop("distops must be a list")
  #}
  if(! (is.null(distops) || is.list(distops))){stop("distops should be NULL or a list")}
  if(! (is.null(fops) || is.list(fops))){stop("fops should be NULL or a list")}

  if (out==TRUE){
    res_out <- perm.test.out(x,y,distops,f,fops,num.perm,diag,exact)
    return(res_out)
  }
  
  
  # Handling function inputs for y
  if (is.function(y)) y <- as.character(substitute(y))
  if (is.character(y)) y <- chartoli(y)
  
  # Calculating observed test statistic
  # One Sample
  if (is.list(y)) {
    if (length(fops) == 0) fops <- NULL
    if (length(distops) == 0) distops <- NULL
    ts.obs <- do.call(f, c(list(x), list(names(y)), distops, fops))
  }
  # Two sample
  if (is.numeric(y)){
    if (length(fops)== 0) fops <- NULL
    ts.obs <- do.call(f, c(list(x,y), fops))
  }
  
  lenx <- length(x)
  ts.random <- vector(mode = "numeric", length = num.perm)
  # Two sample
  if (is.numeric(y)) {
    z <- c(x, y)
    lenz <- length(z)
    if (lenz < 11 & exact == TRUE) {
      require(gtools)
      
      all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, 
                               set = FALSE)
      all.permx <- all.perm[, 1:lenx]
      all.permy <- all.perm[, (lenx + 1):lenz]
      exact.perm <- dim(all.perm)[1]
      for (i in 1:exact.perm) {
        ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
      }
      p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
      c.val <- quantile(ts.random, probs = 0.95)
    } else {
      for (i in 1:num.perm) {
        z1 <- sample(z, size = lenz, replace = FALSE)
        a <- z1[1:lenx]
        b <- z1[(lenx + 1):lenz]
        ts.random[i] <- do.call(f, c(list(a,b), fops))
      }
      p.val <- sum(abs(ts.random) >= abs(ts.obs)) / num.perm
      c.val <- quantile(ts.random, probs = 0.95)
    }
    if (do.plot == TRUE){
      hplot <- hist(ts.random, prob=TRUE)
      hplot
      segments(ts.obs, 0, x1=ts.obs, y1=max(hplot$density))
      segments(c.val, 0, x1=c.val, y1=max(hplot$density), col="red")
    }
    if (diag == TRUE) {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs, 
                  "ts.dist" = ts.random))
    } else {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs))
    }
  }
  # One Sample
  if (is.list(y)) {
    ry <- dist.conv(funname = names(y), type = "r")
    fy <- get(names(y), mode = "function", envir = parent.frame())
    for (i in 1:num.perm) {
      z <- do.call(ry, c(list(lenx), distops))  #(lenx,...)
      ts.random[i] <- do.call(f, c(list(z), list(names(y)), distops, fops))
    }
    p.val <- sum(abs(ts.random) >= abs(ts.obs)) / num.perm
    c.val <- quantile(ts.random, probs = 0.95)
    if (do.plot == TRUE){
      hplot <- hist(ts.random, prob=TRUE, main="Histogram of Permuted Test Statistic (red=Observed TS)")
      hplot
      segments(ts.obs,0,x1=ts.obs,y1=max(hplot$density), col="red")
    }
    if (diag == TRUE) {
      # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
      # statistic
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs, 
                  "ts.dist" = ts.random))
    } else {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs))
    }
  }
}
power.res.onesamp <- function(x, y, distops = NULL, f, fops = NULL, g = perm.test, ...) {
    # Args: 
    # x: numeric matrix 
    # y: either: function name, character naming a function, or
    # list with the format list(qnorm=qnorm)
    dim.x <- dim(x)
    fun <- f
    if (is.character(f)) 
        fun <- get(fun, mode = "function", envir = parent.frame())
    if (is.function(y)) 
        y <- as.character(substitute(y))
    if (is.character(y)) 
        y <- chartoli(y)
    pv <- vector(mode = "numeric", length = dim.x[1])
    for (i in 1:dim.x[1]) {
        a <- g(x[i, ], y, distops, f = fun, fops)
        pv[i] <- a[[1]]
    }
    return(pv)
} 
power.res.twosamp <- function(x, y, distops = NULL, f, fops = NULL, g = perm.test, num.perm = 1001, ...) {
  # Args: 
  # x: numeric matrix 
  # y: either: function name, character naming a function, or
  # list with the format list(qnorm=qnorm)
  dim.x <- dim(x)
  fun <- f
  if (is.character(f)) 
    fun <- get(fun, mode = "function", envir = parent.frame())

  pv <- vector(mode = "numeric", length = dim.x[1])
  for (i in 1:dim.x[1]) {
    a <- g(x[i, ], y[i, ], distops, f = fun, fops, num.perm=num.perm )
    pv[i] <- a[[1]]
  }
  return(pv)
} 
Make_Inter_CDF <- function(x,y,interp=4){
  z <-c(x,y)
  z <- sort(z)
  lenz <- length(z)
  z1 <- seq(1/(lenz+1),lenz/(lenz+1), length.out=lenz)
  #q1 <- quantile(z, probs = z1, type = interp)
  #q_inter <- approxfun(z1, q1, yleft = min(q1), yright = max(q1))
  cdf_inter <- approxfun(z,z1, yleft=0, yright=1)
  return(cdf_inter)
}
Bi_Var_PIT <- function(x,y){
  #com_QFun <- Make_Inter_QFun(x,y)
  inter_CDF <- Make_Inter_CDF(x,y)
  U_x <- inter_CDF(x)
  U_y <- inter_CDF(y)
  return(list(U_x,U_y))
}
Bi_Var_PIT_ks <- function(x,y, alpha=.05){
  # Calculates Bivariate PIT using KS
  a <- Bi_Var_PIT(x,y)
  pval_x <- ks.test(a[[1]],punif,0,1)$p.value
  pval_y <- ks.test(a[[2]],punif,0,1)$p.value
  FWER <- alpha/2
  pvals <- c(pval_x,pval_y)
  if(pval_x < FWER & pval_y < FWER){
    reject <- TRUE
  }else{
    reject <- FALSE
  } 
  names(reject) <- "REJECT NULL?"
  return(list("Reject Null?"=reject, "Pvals"=c(pval_x,pval_y)))
}

scale_two_sample <- function(x,y){
  max_com <- max(x,y)
  min_com <- min(x,y)
  range_x <- max(x) - min(y)
  range_y <- max(y) - min(y)
  range_com <- max_com - min_com
  if(range_x > range_y){
    scaled_x <- (x - min_com) / range_com
    scaled_y <- (y - min_com) / range_com
  } else{
    scaled_x <- (x - min_com) / range_com
    scaled_y <- (y - min_com) / range_com
  }
#   if(range_x > range_y){
#     scaled_x <- (x - min(x)) / range_x
#     scaled_y <- (y - min(x)) / range_x
#   } else{
#     scaled_x <- (x - min(y)) / range_y
#     scaled_y <- (y - min(y)) / range_y
#   }
  return(list(scaled_x, scaled_y))
}

quad.area <- function(x1, x2, y1, y2) {
  t1 <- tri.area(x1, x2, y1)
  t2 <- tri.area(x2, y1, y2)
  area <- t1 + t2
  return(area)
}
tri.area <- function(x, y, z) {
  area <- 0.5 * abs((x[1] - z[1]) * (y[2] - x[2]) - (x[1] - y[1]) * (z[2] - x[2]))
  return(area)
}

#Convenience Function
saveresults <- function(name){
  # name is a string
  temptime <- Sys.time()
  temptime <- gsub(":","", temptime)
  save.image(paste(name, temptime, ".RData", sep=""))
}
