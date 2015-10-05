# Uses Pantelis' cpp code.
# 2 dimensional local weighted least squares smoother. Only local linear smoother is implemented (no higher order, no derivative). 
# bw: bandwidth, a scalar or a vector of length 2.
# kern: kernel used: 'gauss', 'rect', 'gausvar', 'epan' (default), 'quar'
# xin: an n by 2 dataframe or matrix of x-coordinate.
# yin: a vector of y-coordinate.
# win: a vector of weights on the observations. The number of count as in (maybe) raw covariance should be integrated into win.
# xout1: a p1-vector of first output coordinate grid.
# xout2: a p2-vector of second output coordinate grid. If both xout1 and xout2 are unspecified then the output gridpoints are set to the input gridpoints.
# xout: alternative to xout1 and xout2. A matrix of p by 2 specifying the output points.

# Returns a p1 by p2 matrix of fitted values if xout is not specified. Otherwise a vector of length p. 

lwls2dV2 <- function(bw, kern='epan', xin, yin, win=NULL, xout1=NULL, xout2=NULL, xout=NULL, subset=NULL, crosscov = FALSE) {
  if (length(bw) == 1)
    bw <- c(bw, bw)
    
  xin <- matrix(xin, ncol=2)
    
  if (is.null(win))
    win <- rep(1, nrow(xin))
    
  if (!is.null(subset)) {
    xin <- xin[subset, ]
    yin <- yin[subset]
    win <- win[subset]
  }

  if (!is.null(xout1) && !is.null(xout2) && !is.null(xout)) {
    stop('Either xout1/xout2 or xout should be specified.')
  }
  
  if (is.null(xout1)) 
    xout1 <- sort(unique(xin[, 1]))
    
  if (is.null(xout2)) 
    xout2 <- sort(unique(xin[, 2]))

# For passing numerics into the cpp smoother.
  storage.mode(bw) <- 'numeric'
  storage.mode(xin) <- 'numeric'
  storage.mode(yin) <- 'numeric'
  storage.mode(win) <- 'numeric'
  storage.mode(xout1) <- 'numeric'
  storage.mode(xout2) <- 'numeric'

  if (!crosscov){
    ret <- Rmullwlsk(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  } else {
    ret <- RmullwlskCC(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  }
  
  if (!is.null(xout)) {
    ret <- interp2lin(xout1, xout2, ret, xout[, 1], xout[, 2])
  }
  
  return(ret)
}
