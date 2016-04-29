# #' Two dimensional local linear kernel smoother.
# #'
# #' Two dimensional local weighted least squares smoother. Only local linear smoother for estimating the original curve is available (no higher order, no derivative). 
# #' @param bw A scalar or a vector of length 2 specifying the bandwidth.
# #' @param kern Kernel used: 'gauss', 'rect', 'gausvar', 'epan' (default), 'quar'.
# #' @param xin An n by 2 dataframe or matrix of x-coordinate.
# #' @param yin A vector of y-coordinate.
# #' @param win A vector of weights on the observations. 
# #' @param xout1 a p1-vector of first output coordinate grid. Defaults to the input gridpoints if left unspecified.
# #' @param xout2 a p2-vector of second output coordinate grid. Defaults to the input gridpoints if left unspecified.
# #' @param xout alternative to xout1 and xout2. A matrix of p by 2 specifying the output points (may be inefficient if the size of \code{xout} is small).
# #' @param crosscov using function for cross-covariance estimation (Default: FALSE)
# #' @param subset  a vector with the indeces of x-/y-/w-in to be used (Default: NULL)
# #' @return a p1 by p2 matrix of fitted values if xout is not specified. Otherwise a vector of length p corresponding to the rows of xout. 
# #' @export

Lwls2D2 <- function(bw, kern='epan', xin, yin, win=NULL, xout1=NULL, xout2=NULL, xout=NULL, subset=NULL, method = 'SearchTree') {

  # only support epan kernel now.
  stopifnot(kern == 'epan')
  
  if (length(bw) == 1){
    bw <- c(bw, bw)
  }
  if (!is.matrix(xin) ||  (dim(xin)[2] != 2) ){
    stop('xin needs to be a n by 2 matrix')
  }
  # xin <- matrix(xin, ncol=2) # This causes unexcepted/wrong results.
  
  if (is.null(win)){
    win <- rep(1, nrow(xin))
  }  
  if (!is.null(subset)) {
    xin <- xin[subset, ]
    yin <- yin[subset]
    win <- win[subset]
  }
  
  if (!is.null(xout1) && !is.null(xout2) && !is.null(xout)) {
    stop('Either xout1/xout2 or xout should be specified, but not both.')
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
  if (!is.null(xout))
    storage.mode(xout) <- 'numeric' 
  
  if (method == 'SearchTree') {
  # browser()
    xtree <- SearchTrees::createTree(xin)
    ret <- apply(expand.grid(xout1, xout2), 1, function(coord) {
    browser()
      ind <- SearchTrees::rectLookup(xtree, ptOne=coord - bw, ptTwo=coord+bw)
      x <- xin[ind, , drop=FALSE]
      y <- yin[ind]
      w <- win[ind]
      
      llx <- cbind((x[, 1, drop=FALSE] - coord[1]) / bw[1], 
                   (x[, 2, drop=FALSE] - coord[2]) / bw[2])
      temp <- (1 - llx[, 1]^2) * (1 - llx[, 2]^2) * w # simplified; 9/16 is gone
      X <- cbind(1, llx) # simplified; works for original function estimate
      beta <- qr.solve(t(X) %*% diag(temp) %*% X, t(X) %*% diag(temp) %*% y)
      beta[1]
    })
    ret <- matrix(ret, length(xout1), length(xout2))
  } else if (method == 'sort1') {
    ord <- order(xin[, 1])
    xin <- xin[ord, ]
    yin <- yin[ord]
    win <- win[ord]
    # browser()
    ret <- RmullwlskCCsort(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  } else if (method == 'sort2') {
    ord <- order(xin[, 1])
    xin <- xin[ord, ]
    yin <- yin[ord]
    win <- win[ord]
    # browser()
    ret <- RmullwlskCCsort2(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  } else if (method == 'tree') {
    xin[, 1] <- xin[, 1] / bw[1]
    xin[, 2] <- xin[, 2] / bw[2]
    xout1 <- xout1 / bw[1]
    xout2 <- xout2 / bw[2]

    ret <- RmullwlskCCtree(kern, t(xin), yin, win, xout1, xout2)
  } else if (method == 'plain') {
    ret <- RmullwlskCC(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  }
  
  if (!is.null(xout)) {
    ret <- interp2lin(xout1, xout2, ret, xout[, 1], xout[, 2])
  }
  
  return(ret)
}
