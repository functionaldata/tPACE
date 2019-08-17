#' Normalise sparse multivariate functional data
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
#' @examples
#' set.seed(1)
#' n <- 100
#' M <- 51
#' pts <- seq(0, 1, length.out=M)
#' mu <- rep(0, length(pts))
#' sampDense <- MakeGPFunctionalData(n, M, mu, K=1, basisType='sin', sigma=0.01)
#' samp4 <- MakeFPCAInputs(tVec=sampDense$pts, yVec=sampDense$Yn)
#' res4E <- FPCA(samp4$Ly, samp4$Lt, list(error=TRUE))
#' sampN <- GetNormalisedSample(res4E, errorSigma=TRUE)
#'
#' CreatePathPlot(subset=1:20, inputData=samp4, obsOnly=TRUE, showObs=FALSE)
#' CreatePathPlot(subset=1:20, inputData=sampN, obsOnly=TRUE, showObs=FALSE)
#' @export
GetNormalisedSample<- function(fpcaObj, errorSigma = FALSE){
  if (any( 0>=diag(fpcaObj$fittedCov)) ){
    stop("The fitted autocovariance functions appears to have negative or zero values.")
  }

  if (errorSigma){
    sigmaE = fpcaObj$sigma2
  } else {
    sigmaE = 0
  }

  ynorm = mapply(FUN = function(vy, vt){ 
    return( ( vy - approx(y = fpcaObj$mu, x =fpcaObj$workGrid, xout = vt)$y) /
              approx(y = sqrt(sigmaE + diag(fpcaObj$fittedCov)), x =fpcaObj$workGrid, xout = vt)$y)
  }, vy = fpcaObj$inputData$Ly, vt = fpcaObj$inputData$Lt, SIMPLIFY = FALSE)
  return(list(Ly = ynorm, Lt = fpcaObj$inputData$Lt ))
}

#' \code{GetNormalizedSample} is an alias of \code{GetNormalizedSample}
#' @param ... Passed into GetNormalisedSample
#' @export
#' @rdname GetNormalisedSample 
GetNormalizedSample <- function(...) {
  GetNormalisedSample(...)
}
