# Testing Uniformity
A_Uniform <- function(z, k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  F_z <- 1-(1-z)^k
  F_z
}
rA_Uniform <- function(n, k=1){
  if(k<0 ) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z = 1-(1-y)^(1/k)
}


B_Uniform <- function(z,k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  if(length(k) > 1) stop("Haven't vectorized appropriately")
  if(length(z) > 1){
    F_z <- vector(mode="numeric", length=length(z))
    for(i in seq_along(z)){
      if(z[i] <= .5){
        F_z[i] <- 2^(k-1)*z[i]^k
      }
      if(z[i] > .5){
        F_z[i] <- 1 - 2^(k - 1) * (1 - z[i])^k
      }
    }
    return(F_z)
  }
  
  if(z<=.5){
    F_z <- 2^(k - 1) * z^k
    return(F_z)
  }
  if(z >.5){
    F_z <- 1 - 2^(k - 1) * (1 - z)^k
    return(F_z)
  }
}

rB_Uniform <- function(n, k=1){
  if(k < 0) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z <- vector(mode="numeric", length=n)
  for(i in seq_along(z)){
    if(y[i] <= .5){
      z[i] <- (y[i]/(2^(k-1)))^(1/k)
    }
    if(y[i] > .5){
      z[i] <- 1 - ( (1-y[i]) / (2^(k - 1)) )^(1/k)
    }
  }
  return(z)
}




rC_Uniform <- function(n, k=1){
  if(k<0 ) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z <- vector(mode="numeric", length=n)
  for(i in seq_along(z)){
    if(y[i] <= .5){
      z[i] <- .5 - ((.5-y[i])/(2^(k-1)))^(1/k)
    }
    if(y[i] > .5){
      z[i] <- .5 + ((y[i]-.5) / (2^(k-1)))^(1/k)
    }
  }
  return(z)
}







C_Uniform <- function(z,k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  if(length(k)>1) stop("Haven't vectorized appropriately")
  if(length(z) > 1){
    F_z <- vector(mode="numeric", length=length(z))
    for(i in seq_along(z)){
      if( z[i] <=.5)
        F_z[i] <- .5 - 2^(k - 1)*(.5 - z[i])^k
      if(z[i] > .5){
        F_z[i] <- .5 + 2^(k - 1)*(z[i] - .5)^k
      }
    }
    return(F_z)
  }
  if(z <= .5){
    F_z <- .5 - 2^(k-1)*(.5 - z)^k
    return(F_z)
  }
  if(z > .5){
    F_z <- .5 + 2^(k-1)*(z - .5)^k
    return(F_z)
  }
}

# Saw Density
# Want two functions: psaw and rsaw

p2uni <- function(x){
  if(length(x) > 1){
    F_x <- vector(mode="numeric", length=length(z))
    for(i in seq_along(x)){
      if(x[i]<0){
        F_x[i] <- 0
      }
      if(x[i]<1){
        F_x[i] <- (x[i]^2)/2
      }
      if(x[i]>=1 & x[i]<=2){
        F_x[i] <- 2*x[i] - (x[i]^2)/2 - 1
      }
      if(x[i]>2){
        F_x[i] <- 1
      }
    }
    return(F_x)
  }
  
  if(x<0){
    p <- 0
    return(p)
  }
  if(x<1){
    p <- (x^2)/2
    return(p)
  }
  if(x>=1 & x<=2){
    p <- 2*x - (x^2)/2 - 1
    return(p)
  }
  if(x>2){
    p <- 1
    return(p)
  }
  stop("Out of bounds")
}

# Creating Saw Integral
psaw <- function(x){
  if(length(x) > 1){
    F_x <- vector(mode="numeric", length=length(z))
    for(i in seq_along(x)){
      if(x[i] < .5){
        p1 <- 0 
        p2 <- p2uni(4*x[i] - 0)
      }else if(x[i] < 1){
        p1 <- 1
        p2 <- p2uni(4*x[i] - 4*.5)
      }else if(x[i] < 1.5){
        p1 <- 2
        p2 <- p2uni(4*x[i] - 4*1)
      }else{
        p1 <- 3
        p2 <- p2uni(4*x[i] - 4*1.5)
      }
      p <- (p1 + p2)/4
      F_x[i] <- p
    }
    return(F_x)
  }
  if(x < .5){
    p1 <- 0 
    p2 <- p2uni(4*x - 0)
  }else if(x < 1){
    p1 <- 1
    p2 <- p2uni(4*x - 4*.5)
  }else if(x < 1.5){
    p1 <- 2
    p2 <- p2uni(4*x - 4*1)
  }else{
    p1 <- 3
    p2 <- p2uni(4*x - 4*1.5)
  }
  p <- (p1 + p2)/4
  p
}


rsaw <- function(n){
  y <- runif(n)
  z <- vector(mode="numeric", length=n)
  num <- 4
  for(i in seq_along(z)){
    if(y[i] < .25){
      y2 <- y[i] * 4
      if(y2 < .5){
        z[i] <- sqrt(2*y2)
      }else{
        z[i] <- 2-sqrt(2)*sqrt(1-y2)
      }
      z[i] <- 0 + z[i]/4
    }else if(y[i] < .5){
      y2 <- y[i]*4 - 1
      if(y2 < .5){
        z[i] <- sqrt(2*y2)
      }else{
        z[i] <- 2-sqrt(2)*sqrt(1-y2)
      }
      z[i] <- .5 + z[i]/4
    }else if(y[i] < .75){
      y2 <- y[i]*4 - 2
      if(y2 < .5){
        z[i] <- sqrt(2*y2)
      }else{
        z[i] <- 2-sqrt(2)*sqrt(1-y2)
      }
      z[i] <- 1 + z[i]/4
    }else{
      y2 <- y[i]*4 - 3
      if(y2 < .5){
        z[i] <- sqrt(2*y2)
      }else{
        z[i] <- 2-sqrt(2)*sqrt(1-y2)
      }
      z[i] <- 1.5 + z[i]/4
    }
  }
  z
}

# Uniform 
# Uniform Point Mass

rpm <- function(n, p){
  x <- runif(n)
  x2 <- runif(n)
  y <- numeric(n)
  for(i in seq_along(x)){
    if(x[i] < p){
      y[i] <- 0
    }else{
      y[i] <- x2[i]
    }
  }
  y
}

ppm <- function(x, p){
  y <- numeric(length(x))
  for(i in seq_along(x)){
    if(x[i] < 0){
      y[i] <- 0
    }else if(x[i] == 0){
      y[i] <- p
    }else if(x[i]<1){
      y[i] <- p + x[i]*(1-p)
    }else{
      y[i] <- 1
    }
  }
  y
}
