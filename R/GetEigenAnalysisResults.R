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
    
  if (noeig > length(d)) {
    noeig <- length(d);
    warning(sprintf("At most %d number of PC can be selected!", noeig))
  }

  eigenV <- eigenV[, 1:noeig];
  d <- d[1:noeig];

# normalization
  phi <- apply(eigenV, 2, function(x) {
                    x <- x / sqrt(trapz(regGrid, x^2)) 
                    if (x[1] <= x[2])
                      return(x)
                    else
                      return(-x)
  })
  lambda <- gridSize * d;

  fittedCov <- phi %*% diag(lambda) %*% t(phi)
  # Garbage Collection
  gc()
  return(list(lambda = lambda, phi = phi, FVE=FVEObj$FVE[length(FVEObj$FVE)],
              kChoosen=FVEObj$no_opt, fittedCov=fittedCov))
}
