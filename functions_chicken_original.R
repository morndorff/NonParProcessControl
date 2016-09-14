wave.den <- function (x, y, n=2^5, doplot=F) 
{
  # Get cdfs:
  F.x <- ecdf(x)
  F.y <- ecdf(y)
  
  z <- seq(range(x, y)[1], range(x, y)[2], length=n)
  
  F.dwt <- dwt(F.x(z) - F.y(z), wf="haar", n.levels=log(n, 2))
  oc <- unlist(F.dwt)
  oc <- max(abs(oc))
  ks <- max(abs(F.x(z) - F.y(z)))
  
  if(doplot)
  {
    plot(z, F.x(z), type="l", 
         ylim=range(1.1, F.x(z), F.y(z), F.x(z) - F.y(z)))
    lines(z, F.y(z), col=3)
    lines(z, F.x(z) - F.y(z), col=4)
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    segments(z[3], 0, z[3], ks, lwd=2)
    title(paste("oc = ", round(oc, 2), ", ks = ", round(ks, 2), 
                sep=""))
    
  }
  c(oc, ks)
}

wave.den.perm <- function (x=rnorm(10), y=rnorm(10), n=2^5, p=100, doplot=F, doplot1=F) 
{
  z <- c(x, y)
  n.x <- length(x)
  n.z <- length(z)
  b <- numeric(0)
  
  for(i in 1:p)
  {
    perm <- sample(z, n.z, replace=F)
    #		a <- wave.den(perm[1:n.x], perm[(n.x + 1):n.z], n=n, 
    #			doplot=doplot) 
    a <- wave.energy(perm[1:n.x], perm[(n.x + 1):n.z], n=n, 
                     doplot=doplot) 
    b <- rbind(b, a)
  }
  colnames(b) <- c("oc", "ks")
  rownames(b) <- rep("", p)
  
  d.oc <- density(b[, 1])
  d.ks <- density(b[, 2])
  
  pp <- ceiling(0.95 * p)
  oc.crit <- sort(b[, 1])[pp]
  ks.crit <- sort(b[, 2])[pp]
  
  a <- wave.den(x, y, n=n, doplot=doplot)
  
  if(doplot1)
  {
    plot(d.oc, xlim=range(d.oc$x, d.ks$x), 
         ylim=range(d.oc$y, d.ks$y),main="")
    lines(d.ks, col=3)
    abline(v=oc.crit, lty=3)
    abline(v=ks.crit, col=3, lty=3)
    abline(v=a[1])
    abline(v=a[2], col=3)
  }
  
  oc.dec <- ifelse(a[1] < oc.crit, 0, 1)
  ks.dec <- ifelse(a[2] < ks.crit, 0, 1)
  c(oc.dec, ks.dec)
}

wave.den.power <- function (n=2^5, p=100, doplot=F, m=100, doplot1=F) 
{
  b <- numeric(0)
  for(i in 1:m)
  {
    x <- rnorm(100)
    #		y <- rt(100, df=20)
    y <- rnorm(100, sd=1)
    a <- wave.den.perm(x=x, y=y, n=n, p=p, doplot=doplot, 
                       doplot1=doplot1)
    b <- rbind(b, a)
  }
  c(oc=mean(b[, 1]), ks=mean(b[, 2]), frac=mean(b[, 1]) / mean(b[, 2]))
}

wave.energy <- function (x, y, n=2^5, doplot=F) 
{
  # Get cdfs:
  F.x <- ecdf(x)
  F.y <- ecdf(y)
  
  z <- seq(range(x, y)[1], range(x, y)[2], length=n)
  
  F.x.dwt <- dwt(F.x(z), wf="haar", n.levels=log(n, 2))
  F.y.dwt <- dwt(F.y(z), wf="haar", n.levels=log(n, 2))
  x.dwt <- unlist(F.x.dwt)
  y.dwt <- unlist(F.y.dwt)
  #	x.dwt <- (rev(sort(abs(x.dwt))))^2
  #	y.dwt <- (rev(sort(abs(y.dwt))))^2
  
  x.dwt <- rev(sort(abs(x.dwt)))
  y.dwt <- rev(sort(abs(y.dwt)))
  
  xx <- cumsum(x.dwt)
  yy <- cumsum(y.dwt)
  xx <- xx / max(xx)
  yy <- yy / max(yy)
  
  if(doplot){
    plot(xx, type="l", ylim=range(xx, yy, xx - yy))
    lines(yy, col=2)
    lines(xx - yy, col=3)
    abline(h=0, lty=3)
  }
  
  c(sum(abs(xx-yy)), max(abs(xx-yy)))
}