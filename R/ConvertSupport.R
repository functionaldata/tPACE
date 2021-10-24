#' Convert support of a mu/phi/cov etc. to and from obsGrid and workGrid
#' 
#' Convert the support of a given function 1-D or 2-D function from \code{fromGrid} to \code{toGrid}.
#' Both grids need to be sorted. This is an interpolation/convenience function.
#' 
#' @param fromGrid vector of points with input grid to interpolate from
#' @param toGrid vector of points with the target grid to interpolate on
#' @param mu any vector of function to be interpolated
#' @param phi any matrix, each column containing a function to be interpolated 
#' @param Cov a square matrix supported on fromGrid * fromGrid, to be interpolated to toGrid * toGrid.
#' @param isCrossCov logical, indicating whether the input covariance is a cross-covariance. If so then the output is not made symmetric.
#'
#' @export


ConvertSupport <- function(fromGrid, toGrid, mu=NULL, Cov=NULL, phi=NULL, isCrossCov=FALSE) {

  # In case the range of toGrid is larger than fromGrid due to numeric error
  buff <- .Machine$double.eps * max(abs(fromGrid)) * 3
  if (abs(toGrid[1] - fromGrid[1]) < buff)
    toGrid[1] <- fromGrid[1]
  if (abs(toGrid[length(toGrid)] - fromGrid[length(fromGrid)]) < buff)
    toGrid[length(toGrid)] <- fromGrid[length(fromGrid)]
  if ( ( fromGrid[1] - buff  >  toGrid[1]) || 
       ( fromGrid[length(fromGrid)] + buff < toGrid[length(toGrid)]) ) {
    stop("Insufficient size of 'fromGrid'.")}

  if (!is.null(mu)) {# convert mu
    return(MapX1D(fromGrid, mu, toGrid))
  } else if (!is.null(Cov)) {
    mode(fromGrid) <- 'numeric'
    mode(toGrid) <- 'numeric'
    mode(Cov) <- 'numeric'
    gd <- expand.grid(X=toGrid, Y=toGrid)
    ret <- matrix(interp2lin(fromGrid, fromGrid, Cov, gd$X, gd$Y), nrow=length(toGrid))
    if (!isCrossCov) { # ensure that covariance is symmetric
      ret <- 0.5 * (ret + t(ret))
    }
    return(ret)
  } else if (!is.null(phi)) {
    return(MapX1D(fromGrid, phi, toGrid))
  }

}

