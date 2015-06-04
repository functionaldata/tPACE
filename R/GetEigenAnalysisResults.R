# noeig: \approx maxK. The maximal number of eigenfunctions to use for implmentation.
# Values:
# phi: a nRegGrid * no_FVE

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
    warning(sprintf("Warning: at most %d number of PC can be selected!", 
            noeig))
  }

  eigenV <- eigenV[, 1:noeig];
  d <- d[1:noeig];

  phi <- eigenV / sqrt(gridSize);
  lambda <- gridSize * d;

  # # normalized eigen functions
  # eigenV <- apply(eigenV,2, function(x) x / sqrt(trapz(regGrid,(x)^2)) )
  # eigenV <- apply(eigenV,2, function(x) if (x[1] <= x[2]){ x }else{ -x} )
  
  # # interpolate from the normalized eigenfunctions
  # phi <- apply(eigenV,2, function(x) interp1(x= regGrid, y= x, xi=obsGrid, method="spline"))
  # if (noeig ==1){
    # phi <- t(phi)
  # }
  
  # # normalized smoothed eigenfunctions
  # phi  <- apply(phi,2, function(x) x / sqrt(trapz(obsGrid,(x)^2)) )

  return(list(lambda = lambda, phi = phi, FVE=FVEObj$FVE[length(FVEObj$FVE)],
              kChoosen=FVEObj$no_opt))
}
