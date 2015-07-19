fitted.FPCA <-  function (object, ...) {
  # Combine the zero-meaned fitted values (ZMFV) and the interpolated mean (IM)
  # to get the final estimates
  ZMFV = object$xiEst %*% t(object$phi);   
  IM = approx(x= object$obsGrid, y=object$mu, object$workGrid)$y 
  return( t(apply( ZMFV, 1, function(x) x + IM))) 
}
