# TODO: Roxygen documentation

fitted.FPCA <-  function (object, objectDer = NULL, ...) {
  # Combine the zero-meaned fitted values (ZMFV) and the interpolated mean (IM)
  # to get the final estimates
  
  fpcaObj <- object;
  fpcaObjDer <- objectDer;

  if (class(fpcaObj) != 'FPCA'){
    stop("fitted.FPCA() requires an FPCA class object as basic input")
  }

  if ( is.null(objectDer) ){ 
    objToUse <-fpcaObj 
  }  else {
    if(class(objectDer) != "FPCAder"){
      stop("'fitted.FPCA()' has not receive as input a valid FPCAder object.")  
    }
    objToUse <- fpcaObjDer
  }

  ZMFV = fpcaObj$xiEst %*% t(objToUse$phi);   
  IM = approx(x= objToUse$obsGrid, y=objToUse$mu, fpcaObj$workGrid)$y 
  return( t(apply( ZMFV, 1, function(x) x + IM))) 
}
