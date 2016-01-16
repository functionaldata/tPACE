#' Fitted functional sample from FPCA (or FPCAder) object
#' 
#' Combine the zero-meaned fitted values and the interpolated mean to get the fitted values for the trajectories or the first derivatives of these trajectories.
#' 
#' @param object A object of class FPCA returned by the function FPCA().   
#' @param k The integer number of the first k components used for the representation. (default: length(fpcaObj$lambda ))
#' @param der A logical specifying if derivatives should be returned or not (default: FALSE)  
#' @param method The method used to produce the sample of derivatives ('EIG' (default) or 'QUO'). See Liu and Mueller (2009) for more details
#' @param ... Additional arguments
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' fittedY <- fitted(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


fitted.FPCA <-  function (object, k = NULL, der = FALSE, method = NULL, ...) {
  
  fpcaObj <- object;

  if (class(fpcaObj) != 'FPCA'){
    stop("fitted.FPCA() requires an FPCA class object as basic input")
  }

  if( is.null(k) ){
    k = length( fpcaObj$lambda )
  } else {
    if( ( round(k)>=1) && ( round(k) <= length( fpcaObj$lambda ) ) ){
      k = round(k);
    } else {
      stop("'fitted.FPCA()' is requested to use more components than it currently has available. (or 'k' is smaller than 1)")
    }
  }
  
  if( !der ){  
    ZMFV = fpcaObj$xiEst[,1:k, drop = FALSE] %*% t(fpcaObj$phi[,1:k, drop = FALSE]);   
    IM = approx(x= fpcaObj$obsGrid, y=fpcaObj$mu, fpcaObj$workGrid)$y 
    return( t(apply( ZMFV, 1, function(x) x + IM))) 
  } else { #Derivative is not zero
   
     if( is.null(method) ){
      method = 'EIG'
    }

    if( method =='EIG'){
      phi = apply(fpcaObj$phi, 2, getDerivative, t= fpcaObj$workGrid)
      mu = getDerivative(y = fpcaObj$mu, t = fpcaObj$obsGrid)
  
      #if('smoothEIG' == 'FALSE'){
      #  # Smooth very aggressively using splines / Placeholder code
      #  phi = apply(phi,2, function(x) predict(gam(ft~s(t, k= 9), data=data.frame(t=as.vector(fpcaObj$workGrid),ft=x))))
      #  mu = predict(gam(ft~s(t, k= 9), data=data.frame(t=as.vector(fpcaObj$obsGrid),ft=as.vector(mu))))
      #}

      ZMFV = fpcaObj$xiEst[,1:k, drop = FALSE] %*% t(phi[,1:k, drop = FALSE]);
      IM = approx(x= fpcaObj$obsGrid, y=mu, fpcaObj$workGrid)$y
      return( t(apply( ZMFV, 1, function(x) x + IM)))
    }
    if( method == 'QUO'){
      impSample <- fitted(fpcaObj); # Look ma! I do recursion!
      impSampleDer <- t(apply( impSample,1,getDerivative, fpcaObj$workGrid));
      return(impSampleDer)
    }
    warning('You asked for a derivation scheme that is not implemented.')
    return(NULL)
  }
}

getEnlargedGrid <- function(x){
  N <- length(x)
  return (  c( x[1] - 0.1 * diff(x[1:2]), x, x[N] + 0.1 * diff(x[(N-1):N])) )
}

getDerivative <- function(y,t){  # Consider using the smoother to get the derivatives
  if( length(y) != length(t) ){
    stop("getDerivative y/t lengths are unequal.")
  }
  newt = getEnlargedGrid(t) # This is a trick to get first derivatives everywhere
  newy = Hmisc::approxExtrap(x=t, y=y, xout= newt)$y
  return (numDeriv::grad( stats::splinefun(newt, newy) , x = t ) )
}

