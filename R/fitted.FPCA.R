#' Fitted functional sample from FPCA (or FPCAder) object
#' 
#' Combine the zero-meaned fitted values and the interpolated mean to get the final values.
#' 
#' @param object A object of class FPCA returned by the function FPCA().   
#' @param objectDer A object of class FPCAder returned by the function FPCAder(. 
#' @param ... Additional arguments
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' fittedY <- fitted(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


fitted.FPCA <-  function (object, objectDer = NULL, ...) {
  
  fpcaObj <- object;
  fpcaObjDer <- objectDer;

  if (class(fpcaObj) != 'FPCA'){
    stop("fitted.FPCA() requires an FPCA class object as basic input")
  }

  if ( is.null(objectDer) ){ 
    objToUse <-fpcaObj 
  } else {
    if(class(objectDer) != "FPCAder"){
      stop("'fitted.FPCA()' has not receive as input a valid FPCAder object.")  
    }
    objToUse <- fpcaObjDer
  }

  ZMFV = fpcaObj$xiEst %*% t(objToUse$phi);   
  IM = approx(x= objToUse$obsGrid, y=objToUse$mu, fpcaObj$workGrid)$y 
  return( t(apply( ZMFV, 1, function(x) x + IM))) 
}
