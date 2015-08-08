FPCAder <-  function (fpcaObj, ...) {
  # Use FPCA object to get derivative information object 
   
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



