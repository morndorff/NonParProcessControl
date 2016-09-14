# Plotting Functions For making plots based on quantile statistics, 1 sample
plot.ts.1sam <- function(x, y, ..., funname, lenx) {
    # Getting 'y' and converting it to a distribution function
    dist <- strsplit(funname, "q")[[1]][2]
    pdist <- paste("p", dist, sep = "")
    pdist <- get(pdist, mod = "function", envir = parent.frame())
    # 
    y1 <- seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx)
    xcurve <- seq(y(min(y1), ...), y(max(y1), ...), length = 1000)
    ycurve <- pdist(xcurve, ...)
    plot(x, y1, ylim = c(0, 1), xlim = c(min(c(xcurve, x)), max(c(xcurve, x))))
    lines(xcurve, ycurve)
    segments(x, y1, x1 = y(y1, ...), y1 = y1)
    a <- which.max(abs(x - y(y1, ...)))
    segments(x[a], y1[a], x1 = y(y1[a], ...), y1 = y1[a], col = "blue", lwd = 4)
}

# For making plots based on quantile statistics, 2 sample
plot.ts.2sam <- function(x, y, x1, y1, lenx, leny) {
    if (lenx == leny) {
        par(mfrow = c(2, 1))
        # Making Quantile Graph: Note: Because Lengths are the same, don't need
        # iterpolation. (use type=1,for quantiles) For graphical purposes only
        q1 <- quantile(x, probs = x1, type = 1)
        q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1), method = "constant")
        plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
            xlab = "Probs", ylab = "Data")  #quartile plot for z2
        points(y1, y, col = "red")
        points(x1, q1)
        # Making ECDF Graph
        f1 <- ecdf(x)
        plot(f1, xlim = c(min(x[1], y[1]), max(x[lenx], y[lenx])), ylab = "Probs", xlab = "Data", 
            main = "ECDF and Points")  #plot 1st sample points
        points(y, y1, col = "red")  #put 2nd sample on graph
    } else if (lenx > leny) {
        q1 <- quantile(x, probs = x1, type = 4)
        q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
        z <- max(abs(q_inter(y1) - y))
        f <- approxfun(x, x1, yleft = 0, yright = 1, ties = max)  #ties=max makes sure cdf jumps to .4
        par(mfrow = c(2, 1))
        plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
            xlab = "Probs", ylab = "Data")  #quartile plot for z2
        points(y1, y, col = "red")
        plot(f, ylim = c(0, 1), xlim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 
            1), main = "Linearly Interpolated ECDF", xlab = "Data", ylab = "Probs")  #plot of ECDF 
        # Want the limits to include the points, so above w/ maxes and mins is necessary
        lines(ecdf(x))  #draws in original ECDF
        points(y, y1, col = "red")
    } else {
        # So length y>x
        q2 <- quantile(y, probs = y1, type = 4)
        q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
        z <- max(abs(q_inter(x1) - x))
        f <- approxfun(y, y1, yleft = 0, yright = 1, ties = max)  #ties=max makes sure cdf jumps to .4
        par(mfrow = c(2, 1))
        plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
            xlab = "Probs", ylab = "Data")  #quartile plot for z2
        points(x1, x, col = "red")
        plot(f, ylim = c(0, 1), xlim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 
            1), main = "Linearly Interpolated ECDF", xlab = "Data", ylab = "Probs")  #plot of ECDF 
        # Want the limits to include the points, so above w/ maxes and mins is necessary
        lines(ecdf(y))  #draws in original ECDF
        points(x, x1, col = "red")
    }
}

diag.plot <- function(x, y) {
    sortx <- sort(x)
    sorty <- sort(y)
    difxy <- sortx - sorty
    com.xy <- sort(c(x, y))
    # plot(com.xy) A plot to illustrate some of the problem
    x1 <- seq(1/length(x), 1, length.out = length(x))
    y1 <- seq(1/length(y), 1, length.out = length(y))
    plot(ecdf(com.xy), pch = 20)
    segments(sortx, x1, x1 = sorty, y1 = y1, col = "red")
    # want biggest in blue (This code snippet will be useful later)
    a <- which.max(abs(difxy))
    segments(sortx[a], x1[a], x1 = sorty[a], y1 = y1[a], col = "blue", lwd = 4)
}


# QQ Plot
myqq <- function(x, y) {
    # Rank Data
    sx <- sort(x)
    sy <- sort(y)
    # Dealing with unequal sample sizes via linear interpolation May want to do
    # something other than linear interpolation later might be useful to have seperate
    # function
    lenx <- length(sx)
    leny <- length(sy)
    if (lenx > leny) 
        sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- approx(1L:leny, sy, n = lenx)$y
    plot(sx, sy)
    results <- list(sx, sy)
    return(results)
} 

Plot_Quad_Areas <- function (x, y, lenx, coords, xdens, ydens, areas=NULL) {
  # Plotting the Quad_OS_Area_TS
  # x_density comes from inside the TS function
  x_points <- seq(min(x)-sd(x), max(x)+sd(x), length.out=500)
  x_out <- sapply(x_points, function(x) xdens(x))
  
  y_points <- seq(min(y)-sd(y), max(y)+sd(y), length.out=500)
  y_out <- sapply(y_points, function(x) ydens(x))
  
  plot(x_points,x_out, xlim = c(min(x[1]-sd(x), y[1]-sd(y)), max(x[lenx]+sd(x), y[lenx]+sd(y))), ylab = "Probs", xlab = "Data", 
       main = "Two Estimated ECDFS, Areas in Blue", type="l")  #plot 1st sample points
  lines(y_points,y_out,col="red")
  
  plot_polygons <- function(liCoords, color="blue"){
    plotx <- sapply(liCoords,function(x) x[1])  
    ploty <- sapply(liCoords,function(x) x[2])  
    my_plot <- polygon(plotx,ploty, col=color)
    my_plot
  }
  sapply(coords, plot_polygons)
  if(is.null(areas)==FALSE){
    big_area <- as.numeric(which.max(areas))
    plot_polygons(coords[[big_area]], color="red")
  }
}

#Plot_Quad_Quantiles
