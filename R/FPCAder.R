#' Functional Principal Component Analysis Derivatives
#' 
#' Derivative FPCA for dense or sparse functional data. 
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA(). 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' resder <- FPCAder(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


FPCAder <-  function (fpcaObj, ...) {
  # Use FPCA object to get derivative information object 
   
  if (class(fpcaObj) != 'FPCA'){
    stop("FPCAder() requires an FPCA class object as basic input")
  }
  
  fpcaObjDer <- list( 
    phi = apply(fpcaObj$phi, 2, getDerivative, t= fpcaObj$workGrid),
    mu = getDerivative(y = fpcaObj$mu, t = fpcaObj$obsGrid), 
    obsGrid = fpcaObj$obsGrid)
     
  class(fpcaObjDer) <- 'FPCAder'
  return(fpcaObjDer) 
} 

getEnlargedGrid <- function(x){ 
  N <- length(x)
  return (  c( x[1] - 0.1 * diff(x[1:2]), x, x[N] + 0.1 * diff(x[(N-1):N])) )
}

getDerivative <- function(y,t){ 
  newt = getEnlargedGrid(t)
  newy = approxExtrap(x=t, y=y, xout= newt)$y 
  return (numDeriv::grad( splinefun(newt, newy) , x = t ) )
}



