# 2 dimensional local weighted least squares smoother. Only local linear smoother is implemented (no higher order, no derivative). 
# bw: bandwidth, a scalar.
# kern: kernel used: 'gauss', 'rect', 'tria', 'epan' (default), 'quar', 'trwt', or 'tcub'. (See \url{http://en.wikipedia.org/wiki/Kernel_%28statistics%29#In_non-parametric_statistics})
# xin: an n by 2 dataframe or matrix of x-coordinate.
# yin: a vector of y-coordinate.
# win: a vector of weights on the observations. The number of count as in (maybe) raw covariance should be integrated into win.
# xout1: a p1-vector of first output coordinate grid.
# xout2: a p2-vector of second output coordinate grid. If both xout1 and xout2 are unspecified then the output gridpoints are set to the input gridpoints.
# xout: alternative specification of output points. a matrix of two columns.
# returnFit: If TRUE, return the fitted locfit object for the convenience of future operations such as gcv. Use predict(fit) to obtain the fitted values.
# ...: pass on to locfit
# Returns a p1 by p2 matrix of fitted values.
# Difference from mullwlsk: The local window is elliptical. There is no nwe argument. The kernels supported are different. No local bandwidth choice is supported as it is never actually called in Matlab PACE. 
lwls2d <- function(bw, kern='epan', xin, yin, win=NULL, xout1=NULL, xout2=NULL, xout=NULL, returnFit=FALSE, scale=FALSE, ...) {
  datin <- data.frame(cbind(xin, as.numeric(yin)))
  names(datin) <- c('x1', 'x2', 'y')

  if (missing(xout)) {
    if (missing(xout1) && missing(xout2))
      xout <- as.matrix(datin[, -3])
    else
      xout <- as.matrix(expand.grid(xout1, xout2))
  }
  if (class(xout)[1] != 'lf_evs')
    colnames(xout) <- names(datin)[1:2]

  if (missing(win) || is.null(win))
    win <- rep(1, nrow(datin))

  # decide space allocation
  uniqPts <- sort(unique(datin$x1))
  r <- diff(range(uniqPts))
  # maxk <- round((2 * bw / r * length(uniqPts))^2 * 4)
  maxKLocfit <- 2000
  
  
  if (class(xout)[1] == 'lf_evs') { # use locfit evaluation structure 
    fit <- locfit(y ~ lp(x1, x2, h=bw, deg=1, scale=scale), data=datin, weights=win, kern=kern, ev=xout, maxk=maxKLocfit, ...)
    ret <- fitted(fit, datin)
  } else {
    # browser()
    fit <- locfit(y ~ lp(x1, x2, h=bw, deg=1, scale=scale), data=datin, weights=win, kern=kern, maxk=maxKLocfit, ...)
    # plot(fit)
    ret <- predict(fit, xout)
  }
  
  if (!missing(xout1) && !missing(xout2)) {
    ret <- matrix(ret, length(xout1), length(xout2))
  }

  if (returnFit)
  return(fit)
  else {
    return(ret)
  }
}
