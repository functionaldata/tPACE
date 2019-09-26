#' Two dimensional local linear kernel smoother to target derivatives.
#'
#' Two dimensional local weighted least squares smoother. Only a local linear smoother for estimating the original curve is available (no higher order)
#' @param bw A scalar or a vector of length 2 specifying the bandwidth.
#' @param kern Kernel used: 'gauss', 'rect', 'gausvar', 'epan' (default), 'quar'.
#' @param xin An n by 2 data frame or matrix of x-coordinate.
#' @param yin A vector of y-coordinate.
#' @param win A vector of weights on the observations. 
#' @param xout1 a p1-vector of first output coordinate grid. Defaults to the input gridpoints if left unspecified.
#' @param xout2 a p2-vector of second output coordinate grid. Defaults to the input gridpoints if left unspecified.
#' @param xout alternative to xout1 and xout2. A matrix of p by 2 specifying the output points (may be inefficient if the size of \code{xout} is small).
#' @param npoly The degree of polynomials (include all \eqn{x^a y^b} terms where \eqn{a + b <= npoly})
#' @param nder1 Order of derivative in the first direction
#' @param nder2 Order of derivative in the second direction
#' @param crosscov using function for cross-covariance estimation (Default: TRUE)
#' @param subset  a vector with the indices of x-/y-/w-in to be used (Default: NULL)
#' @param method should one try to sort the values xin and yin before using the lwls smoother? if yes ('sort2' - default for non-Gaussian kernels), if no ('plain' - fully stable; de)
#' @return a p1 by p2 matrix of fitted values if xout is not specified. Otherwise a vector of length p corresponding to the rows of xout. 
#' @export

Lwls2DDeriv <- function(
  bw, kern='epan', xin, yin, win=NULL, xout1=NULL, xout2=NULL, xout=NULL, 
  npoly=1L, nder1=0L, nder2=0L, subset=NULL, 
  crosscov = TRUE, method = 'sort2'
) {

  # only support epan kernel now.
  # stopifnot(kern == 'epan')
  
  if (length(bw) == 1){
    bw <- c(bw, bw)
  }
  if (is.data.frame(xin)) {
    xin <- as.matrix(xin)
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
  storage.mode(npoly) <- 'integer'
  storage.mode(nder1) <- 'integer'
  storage.mode(nder2) <- 'integer'
  
  if (method == 'sort2') {
    ord <- order(xin[, 1])
    xin <- xin[ord, ]
    yin <- yin[ord]
    win <- win[ord]
    # browser()
    ret <- RmullwlskUniversalDeriv(bw, kern, t(xin), yin, win, 
      xout1, xout2, npoly, nder1, nder2, FALSE, !crosscov)
  } else if (method == 'plain') { # MAYBE IMPROVE THIS
    ret <- RmullwlskCC(bw, kern, t(xin), yin, win, xout1, xout2, FALSE)
  }

  if (!is.null(xout)) {
    ret <- interp2lin(xout1, xout2, ret, xout[, 1], xout[, 2])
  }
  
  return(ret)
}
