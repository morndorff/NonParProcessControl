# Testing Functions, October 2, 2014
perm.test2 <- function(x, y, distops = NULL, f, fops = NULL, num.perm = 2001, diag = FALSE, 
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
  if (is.null(distops)==FALSE){
    if(is.list(distops)==FALSE) stop("distops must be a list")
  }
  if (is.null(fops)==FALSE){
    if(is.list(fops)==FALSE) stop("fops must be a list")
  }
  if (out==TRUE){
    res_out <- perm.test.out(x,y,distops,f,fops,num.perm,diag,exact)
    return(res_out)
  }
  lenx <- length(x)
  
  # Handling function inputs for y
  if (is.function(y)) 
    y <- as.character(substitute(y))
  
  # Calculating observed test statistic
  # One Sample
  if (is.character(y)) {
    y <- chartoli(y)
    if (length(fops) == 0) 
      fops <- NULL
    if (length(distops) == 0) 
      distops <- NULL
    ts.obs <- do.call(f, c(list(x), list(names(y)), distops, fops))
  }
  # Two sample
  if (is.numeric(y)){
    if (length(fops)== 0)
      fops <- NULL
    if (is.null(fops[[1]])) # MUST BE EDITED ASAP. ONLY LETS ONE OPTION IN FOPS
      fops <- NULL
    ts.obs <- do.call(f, c(list(x,y), fops))
  }
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
      matPerm <- replicate(num.perm, sample(z, size=lenz, replace=FALSE), simplify="array")
      ts.random <- apply(matPerm, 2, function(x) do.call(f, c(list(x[1:lenx], x[(lenx+1):lenz]), fops)))

#       for (i in 1:num.perm) {
#         z1 <- sample(z, size = lenz, replace = FALSE)
#         a <- z1[1:lenx]
#         b <- z1[(lenx + 1):lenz]
#         ts.random[i] <- do.call(f, c(list(a,b), fops))
#       }
      p.val <- sum(abs(ts.random) >= abs(ts.obs)) / num.perm
      c.val <- quantile(ts.random, probs = 0.95)
    }
    if (do.plot == TRUE){
      hplot <- hist(ts.random, prob=TRUE)
      hplot
      segments(ts.obs, 0, x1=ts.obs, y1=max(hplot$density))
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
      hplot <- hist(ts.random, prob=TRUE, main="Histogram of Permuted Test Statistic (line=Observed TS)")
      hplot
      segments(ts.obs,0,x1=ts.obs,y1=max(hplot$density))
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