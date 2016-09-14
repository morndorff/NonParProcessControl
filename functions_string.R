# String Manipulation Functions
dist.conv <- function(funname, type) {
    # Getting a function name and converting it to a distribution function or random
    # variable genderation 
    # Input 
    # funname: a function name, such as qnorm 
    # type: any character either p (returns pnorm) or r (returns rnorm) 
    # Returns: qnorm -> rnorm, qnorm -> pnorm, etc.
    if (is.function(funname)) {
        funname <- as.character(substitute(funname))
    }
    if (is.character(funname)) {
        dist <- substring(funname, 2)
        tdist <- paste(type, dist, sep = "")
        tdist <- get(tdist, mod = "function", envir = parent.frame())
    }
    return(tdist)
}

dist.conv.str <- function(funname, type) {
    # Getting a function name and converting it to a distribution function or random
    # variable generation 
    # Input funname: a function name, such as qnorm 
    # type: any character either p (returns pnorm) or r (returns rnorm) 
    # Returns: qnorm -> rnorm, qnorm -> pnorm, etc.
    if (is.function(funname)) {
        funname <- as.character(substitute(funname))
    }
    if (is.character(funname)) {
        dist <- substring(funname, 2)
        tdist <- paste(type, dist, sep = "")
    }
    return(tdist)
}

funtochar <- function(x) {
    x <- as.character(substitute(x))
    return(x)
}

make_sample <- function(n, dist, params) {
    # Args: 
    # dist: 'norm', 'unif', etc. 
    # n: sample size 
    # params: list(min=0,max=2) 
    # Ex: sam <- make_sample(n=100,dist=unif,params=list(min=0,max=2)) 
    # sam <- make_sample(n=1000,dist=norm,params=list(mean=0,sd=2))
    dist <- eval.parent(dist[[1]], n = 1)
    dist <- as.character(substitute(dist))
    rdist <- paste("r", dist, sep = "")
    rdist <- get(rdist, mode = "function", envir = parent.frame())
    sample <- do.call(rdist, c(list(n = n), params))
    return(sample)
}

funtoli <- function(fun) {
    a <- list(fun)
    names(a) <- as.character(substitute(fun))
    return(a)
}

chartoli <- function(char) {
    fun <- get(char, mode = "function", envir = parent.frame())
    a <- list(fun)
    names(a) <- char
    return(a)
} 
