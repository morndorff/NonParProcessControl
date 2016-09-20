# Functions for Power Simulation Studies (ks, cramer von mises, etc.)

# K test (Just outputting statistic)
ks.res.simp <- function(x, y, ..., alternative = c("two.sided", "less", "greater"), 
    exact = NULL) {
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1L) 
        stop("not enough 'x' data")
    PVAL <- NULL
    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        if (is.null(exact)) 
            exact <- (n.x * n.y < 10000)
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            } else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), greater = max(z), 
            less = -min(z))
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below that of y", 
            greater = "the CDF of x lies above that of y")
    } else {
        if (is.list(y)) 
            y <- names(y)
        if (is.character(y)) {
            y <- dist.conv.str(y, type = "p")
            y <- get(y, mode = "function", envir = parent.frame())
        }
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if (is.null(exact)) 
            exact <- (n < 100) && !TIES
        x <- y(sort(x), ...) - (0:(n - 1))/n
        STATISTIC <- switch(alternative, two.sided = max(c(x, 1/n - x)), greater = max(1/n - 
            x), less = max(x))
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", greater = "D^+", less = "D^-")
    
    RVAL <- list(statistic = STATISTIC, alternative = nm_alternative, method = METHOD, 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(STATISTIC)
}

ks.two.quick.cont <- function(x,y){
    x <- sort(x)
    y <- sort(y)
    com <- c(x,y)
    ranks <- rank(com)
    lenx <- length(x)
    leny <- length(y)
    i <- 1:lenx
    j <- 1:leny
    d1 <- i/lenx - (ranks[1:lenx] - i)/leny
    d2 <- (ranks[(lenx+1) : (lenx + leny)] - j)/lenx -(j/leny)
    STAT <- max(abs(d1), abs(d2))
    return(STAT)  
}

cvmtwo.res <- function(x,y){
 x <- sort(x)
 y <- sort(y)
 com <- c(x,y)
 ranks <- rank(com)
 lenx <- length(x)
 leny <- length(y)
 i <- 1:lenx
 j <- 1:leny
 U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
 STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
 return(STAT)  
}



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
  #if (is.null(distops)==FALSE){
  #  if(is.list(distops)==FALSE) stop("distops must be a list")
  #}
  library(doMC)
  library(foreach)
  registerDoMC(2)
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
#       for (i in 1:num.perm) {
#         z1 <- sample(z, size = lenz, replace = FALSE)
#         a <- z1[1:lenx]
#         b <- z1[(lenx + 1):lenz]
#         ts.random[i] <- do.call(f, c(list(a,b), fops))
#       }
#       foreach(i=1:num.perm) %do% {
#         z1 <- sample(z, size = lenz, replace = FALSE)
#         a <- z1[1:lenx]
#         b <- z1[(lenx + 1):lenz]
#         ts.random[i] <- do.call(f, c(list(a,b), fops))
#       }
      z_mat <- replicate(num.perm,sample(z,size=lenz, replace=FALSE))
      ts.random <- apply(z_mat,2, function(x) do.call(f, c( list(x[1:lenx], x[(lenx +1):lenz], fops ))))
      # Slower
      
      
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

cvm.res <- function(x, y, ..., test=FALSE){
  # x: numeric vector
  # y: numeric vector or prob distribution e.g. "pnorm" 
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    x <- sort(x)
    y <- sort(y)
    com <- c(x,y)
    ranks <- rank(com)
    lenx <- length(x)
    leny <- length(y)
    i <- 1:lenx
    j <- 1:leny
    U <- lenx * sum((ranks[1:lenx] - i)^2) + leny * sum((ranks[(lenx + 1):(lenx + leny)] - j)^2)
    STAT <- U / ((lenx * leny)*(lenx + leny)) - (4 * leny * lenx - 1)/(6 * (leny + lenx))
    return(STAT)
  }
  # One Sample Test  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  F_x <- y(x, ...)
  i <- 1:lenx
  STAT <- (1/(12*lenx)) + sum( (F_x -(1:lenx - .5)/ lenx)^2)
  if(test){
    # print(F_x)
    # print(F_x - (1:lenx -.5)/lenx)
    # return(F_x - (1:lenx -.5)/lenx)
    STAT <- (1/(12*lenx)) + sum((F_x - (1:lenx - .5)/lenx)^2)
    return(STAT)
  }
  STAT
}

ad.res <- function(x,y, ...){
  # Calculates "A"
  # Two Sample
  if(is.numeric(y)){
    # From 'A two-sample Anderson-Darling rank statistic
    com <- c(x,y)
    lenx <- length(x) #n
    leny <- length(y) #m
    len_com <- lenx + leny
    i <- 1:(len_com-1)
    fx <-ecdf(x)
    M <- lenx * fx(sort(com))
    print(M)
    plot(M)
    M <- M[-len_com]
    STAT <- 1/(lenx*leny)*sum((M * len_com -lenx *i)^2/(i*(len_com-i)))
    return(STAT)
  }
  
  # One Sample Variant First
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  x <- sort(x)
  F_x <- y(x, ...)
  lenx <- length(x)
  ln_F_x <- log(F_x)
  S <- 0
  for(i in 1:lenx){
    S <- (2 * i - 1) * (ln_F_x[i] + log(1 - F_x[lenx + 1 -i])) + S
  }
  S <- S/lenx
  STAT <- -lenx - S
}

ad.check <- function(x,y, ...){
  if(is.numeric(y)){
    # From 'A two-sample Anderson-Darling rank statistic
    x <- sort(x)
    y <- sort(y)
    ranks(c(x,y))
    
    
    com <- c(x,y)
    lenx <- length(x)
    leny <- length(y)
    len_com <- lenx + leny
    i <- 1:(len_com-1)
    fx <-ecdf(x)
    M <- lenx *fx(sort(com))
    M <- M[-len_com]
    STAT <- 1/(lenx*leny)*sum((M * len_com -leny *i)^2/(i*(len_com-i)))
    return(STAT)
  }
  # One Sample Variant First
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  x <- sort(x)
  F_x <- y(x, ...)
  lenx <- length(x)
  Fn <- ecdf(x)
  Fn_x <- Fn(x)
  num <- (Fn_x - F_x)^2
  denom <- (1 - F_x) * (F_x)
  STAT <- sum(num/denom)
  STAT
}

