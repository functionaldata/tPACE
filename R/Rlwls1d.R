#' One dimensional local linear kernel smoother
#' 
#' One dimensional local linear kernel smoother for sparse data.
#'
#' @param bw Scalar holding the bandwidth
#' @param kernel_type Character holding the kernel type (see ?FPCA for supported kernels)
#' @param win Vector of length N with weights
#' @param xin Vector of length N with measurement points
#' @param yin Vector of length N with measurement values
#' @param xout Vector of length M with output measurement points
#' @param npoly Scalar (integer) degree of polynomial fitted (default 1)
#' @param nder Scalar (integer) degree of derivative fitted (default 0)
#'
#' @return Vector of length M with measurement values at the the point speficied by 'xout'
#'
#' @export


Rlwls1d <- function( bw, kernel_type, win, xin, yin, xout, npoly = 1, nder = 0){
  return( CPPlwls1d(bw= as.numeric(bw), kernel_type = kernel_type, npoly= npoly, nder= nder, 
                  xin = as.numeric(xin), yin= as.numeric(yin), xout= as.numeric(xout), win = as.numeric(win)))
}

