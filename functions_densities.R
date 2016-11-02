# Probability Densities

standard_chi <- function(n, df, center_mean=0, center_sd=1){
 omean <- df
 ovar <- 2*df
 
 STAT <- rchisq(n=n, df=df)
 STAT <- (STAT - (omean))/ ((sqrt(ovar)) * 1/center_sd) + center_mean
}

standard_normal <- function(n, mean, sd, center_mean=0, center_sd=1){
  omean <- mean
  ovar <- sd^2
  STAT <- rnorm(n=n, mean=mean, sd=sd)
  STAT <- (STAT - (omean))/ ((sqrt(ovar)) * 1/center_sd) + center_mean
}

standard_t <- function(n, df, center_mean=0, center_sd=1){
  omean <- 0
  if( df<2) stop("I don't support df less than 2!")
  ovar <- df/(df-2)
  STAT <- rt(n=n, df=df)
  STAT <- (STAT - (omean))/ ((sqrt(ovar)) * 1/center_sd) + center_mean
}
