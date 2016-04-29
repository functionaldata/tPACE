#' Normalise sparse functional sample
#'
#' Normalise sparse functional sample given in an FPCA object
#'
#' @param fpcaObj An FPCA object.
#' @param errorSigma Indicator to use sigma^2 error variance when normalising the data (default: FALSE)
#'
#' @return A list containing the normalised sample 'y' at times 't'
#'
#' @references
#' \cite{Chiou, Jeng-Min and Chen, Yu-Ting and Yang, Ya-Fang. "Multivariate Functional Principal Component Analysis: A Normalization Approach" Statistica Sinica 24 (2014): 1571-1596}
#' @export
GetNormalisedSample<- function(fpcaObj, errorSigma = FALSE){
  if (any( 0>=diag(fpcaObj$smoothedCov)) ){
    stop("The autocovariance functions appears to have negative or zero values.")
  }

  if (errorSigma){
    sigmaE = fpcaObj$sigma2
  } else {
    sigmaE = 0
  }

  ynorm = mapply(FUN = function(vy, vt){ 
    return( ( vy - approx(y = fpcaObj$mu, x =fpcaObj$workGrid, xout = vt)$y) /
              approx(y = sqrt(sigmaE + diag(fpcaObj$smoothedCov)), x =fpcaObj$workGrid, xout = vt)$y)
  }, vy = fpcaObj$inputData$Ly, vt = fpcaObj$inputData$Lt, SIMPLIFY = FALSE)
  return(list(Ly = ynorm, Lt = fpcaObj$inputData$Lt ))
}
