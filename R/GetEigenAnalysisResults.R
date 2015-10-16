# noeig: \approx maxK. The maximal number of eigenfunctions to use for implmentation.
# Values:
# phi: a nRegGrid * no_FVE
# The input smoothCov is possibly truncated.

GetEigenAnalysisResults <- function(smoothCov, regGrid, optns) {
#   noeig \approx  maxK 
  noeig <- optns$maxK
  FVEthreshold <- optns$FVEthreshold
  verbose <- optns$verbose
  
  gridSize <- regGrid[2] - regGrid[1]
  numGrids <- nrow(smoothCov)
  
  FVEObj <- no_FVE(smoothCov, FVEthreshold=FVEthreshold, returnEVec=TRUE, verbose=verbose)
  #optns.v0 <- t(seq(0.1,0.9,length.out = numGrids))
  d <- FVEObj$lambda
  eigenV <- FVEObj$eVec
  
  if(!is.matrix(eigenV)){
    eigenV <- matrix(eigenV); # In case it is a vector
  }

  if (noeig > length(d)) {
    noeig <- length(d);
    warning(sprintf("At most %d number of PC can be selected!", noeig)) 
    # the warning shows total number of pos eigenvalues from fitted cov, 
    # not getting this means we only take first maxk components
  }
  # after checking for maxk, get updated eigenvalues and eigenfns
  eigenV <- eigenV[, 1:noeig];
  d <- d[1:noeig];
  # thresholding for corresponding FVE option 
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) * 100 - 0.001 # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVEthreshold*100)) # final number of component chosen based on FVE

  if(!is.matrix(eigenV)){
    eigenV <- matrix(eigenV); # In case it is a vector
  }

# normalization
  phi <- apply(eigenV, 2, function(x) {
                    x <- x / sqrt(trapzRcpp(regGrid, x^2)) 
                    if (x[1] <= x[2])
                      return(x)
                    else
                      return(-x)
  })
  lambda <- gridSize * d;

  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)
  # Garbage Collection
  gc()
  return(list(lambda = lambda[1:no_opt], phi = as.matrix(phi[,1:no_opt]), cumFVE = FVE,
              kChoosen=no_opt, fittedCov=fittedCov))
}
