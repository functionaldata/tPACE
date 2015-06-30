# This function is used to convert the support of mu/phi/cov etc, to and from obsGrid and regGrid.
# If one wants to convert both mu and phi, it should be called one at a time.
# mu: any vector of a function
# phi: any matrix, each column containing a function to be interpolated.
# Cov: any matrix supported on fromGrid * fromGrid, to be interpolated to toGrid * toGrid. 

ConvertSupport <- function(fromGrid, toGrid, mu=NULL, Cov=NULL, phi=NULL) {

  if (!is.null(mu)) {# convert mu
    return(mapX1d(fromGrid, mu, toGrid))
  } else if (!is.null(Cov)) {
    gd <- meshgrid(toGrid)
    ret <- matrix(interp2lin(fromGrid, fromGrid, Cov, gd$X, gd$Y), nrow=length(toGrid))
    return(ret)
  } else if (!is.null(phi)) {
    return(mapX1d(fromGrid, phi, toGrid))
  }

}

