getEigens <- function(userCov,obsGrid,regGrid,noeig, varargin){
#   function [lambda, phi, eigen, noeig] = getEigens(userCov,obsGrid,regGrid,noeig, varargin)
#   noeig \approx  maxK

  h=diff(range(regGrid))/(length(regGrid)-1);

  numGrids = nrow(userCov)
  #optns.v0 = t(seq(0.1,0.9,length.out = numGrids))
  eigsOutput = eigs(userCov, k = numGrids-2, which = "LR")
  d = eigsOutput$values
  eigenV = eigsOutput$vectors
  # at most ngrid-2 eigenvalues can be obtained for nonsymmetric or complex problems
  # "LM" corresponds to Largest Magnitude, another option may be "LR": Largest Real part
  
  idx = which(is.complex(d)) # to remove any imaginary eigenvalues
  if(invalid(idx) == FALSE){ # if there are any imaginary eigenvalues
    stop(sprintf("%d eigenvalues are complex. The estimated auto-covairance surface is not symmetric!",length(idx)))
  }
  
  idx = which(d <= 0) # to remove nonpositive eigenvalues
  if(invalid(idx) == FALSE){ # if there are any nonpositive eigenvalues
    warning(sprintf("Warning: %d real eigenvalues are negative or zero and are removed!",length(idx)))
    d = d[d > 0]
    eigenV = eigenV[,d>0]
  }
    
  if (noeig > length(d)){
    noeig = length(d);
    warning(sprintf("Warning: at most %d number of PC can be selected!", (noeig)  ))  
  }

  eigenV = eigenV[,1:noeig];
  d = d[1:noeig];

  eigenV = eigenV/sqrt(h);
  lambda = h * d;

  # normalized eigen functions
  eigenV = apply(eigenV,2, function(x) x / sqrt(trapzRcpp(regGrid,(x)^2)) )
  eigenV = apply(eigenV,2, function(x) if (x[1] <= x[2]){ x }else{ -x} )
  
  # interpolate from the normalized eigenfunctions
  phi = apply(eigenV,2, function(x) spline(x= regGrid, y= x, xout=obsGrid)$y)
  if (noeig ==1){
    phi = t(phi)
  }
  
  # normalized smoothed eigenfunctions
  phi  = apply(phi,2, function(x) x / sqrt(trapzRcpp(obsGrid,(x)^2)) )

  return( list(lambda = lambda, phi = phi, eigen = eigenV, noeig = noeig) )
}

