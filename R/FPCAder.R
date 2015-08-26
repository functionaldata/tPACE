#' Functional Principal Component Analysis Derivatives
#' 
#' Derivative FPCA for dense or sparse functional data. 
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA(). 
#' @param variant A string specifying the methodology used ('simple', 'QUO'; default: 'simple')
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


FPCAder <-  function (fpcaObj, variant = 'simple') {
  # Use FPCA object to get derivative information object 
   
  if (class(fpcaObj) != 'FPCA'){
    stop("FPCAder() requires an FPCA class object as basic input")
  }
  if( variant == 'simple') {
    fpcaObjDer <- list( 
      phi = apply(fpcaObj$phi, 2, getDerivative, t= fpcaObj$workGrid),
      mu = getDerivative(y = fpcaObj$mu, t = fpcaObj$obsGrid), 
      obsGrid = fpcaObj$obsGrid)
       
    class(fpcaObjDer) <- 'FPCAder'
    return(fpcaObjDer) 
  } else if ( variant == 'QUO' ){
    impSample <- fitted(fpcaObj);
    impSampleDer <- t(apply( impSample,1,getDerivative, fpcaObj$workGrid));
    N = dim(impSample)[1];
    M = dim(impSample)[2];
    L = makePACEinputs(IDs = rep(1:N,each=M), tVec=rep(fpcaObj$workGrid,N), as.vector(t(impSampleDer)))
    prefpcaObjDer = FPCA(y= L$Ly, t= L$Lt)
    fpcaObjDer = list( 
      phi = prefpcaObjDer$phi, mu =  prefpcaObjDer$mu, obsGrid = fpcaObj$obsGrid)     

    class(fpcaObjDer) <- 'FPCAder'
   return(fpcaObjDer)
  } else {
    stop("Invalid FPCAder variant requested.")
    return( NULL )
  }
} 

getEnlargedGrid <- function(x){ 
  N <- length(x)
  return (  c( x[1] - 0.1 * diff(x[1:2]), x, x[N] + 0.1 * diff(x[(N-1):N])) )
}

getDerivative <- function(y,t){
  if( length(y) != length(t) ){
    stop("getDerivative y/t lengths are unequal.")    
  }
  newt = getEnlargedGrid(t)
  newy = Hmisc::approxExtrap(x=t, y=y, xout= newt)$y 
  return (numDeriv::grad( stats::splinefun(newt, newy) , x = t ) )
}



