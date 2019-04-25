# Approximate the minimum bandwidth choice for% for the covariance function. 
# In stead of the getMinb.m functionality we can garantee minimal number of neighboring points in here.
# TODO: distMat is memory inefficient.
GetMinb <- function(t, obsGrid, dataType='Sparse', npoly=1, minUniqPts=3, minPts=6, legacyCode = FALSE) {
  
 
  if( legacyCode ){
    if (dataType == 'Sparse') {
      dstar <- Minb(obsGrid, 2 + npoly) # rough 1D initial value 
      n_obs = length(obsGrid);
      tmp1 =  matrix( rep(0, n_obs^2), ncol = n_obs)
      
      # Find the pair against which we have measurements in the same curve
      for (i in 1:length(t)){
        idx = match( t[[i]], obsGrid)
        tmp1[idx, idx] = 1
      }
      res = tmp1 - diag(n_obs);
      # First and last timepoint are always considered observed
      res[c(1, n_obs),] = 1;  
      ids = matrix(res > 0); 
      b = matrix( rep(obsGrid, n_obs), nrow=n_obs);
      # Use half of the largest difference between two consequative points in the same 
      # as curve as your candidate bandwith. We do no worry about the difference
      # between to [t_j(end) - t_{1+j}(1)] because this will be negative. This bandwidth tends to be conservative (too large).
      # dstar = max(dstar, max(diff(b[ids])/2)); # Original code
      dstar = max(dstar, quantile( diff(b[ids]), 0.95)/2 ); # Fix to avoid outliers
    } else if (dataType == 'RegularWithMV') {
      dstar <- Minb(obsGrid, 1 + npoly) * 2;
    } else if (dataType == 'Dense') {
      dstar = Minb(obsGrid, 2 + npoly) * 1.5;
    }   
    return(dstar)
  }
  
  dstar = BwNN(t, k= 2 + npoly, onlyCov = TRUE)['cov']
  
  return(dstar)
}
