# This function obtains the sample covariance matrix at observed grid
# for dense regular functional data

######
# Input:
######  
#  ymat: n by p matrix of dense regular functional data
#  mu: p-dim vector, estimated cross-sectional mean
#  optns: options for FPCA function
######
# Output: 
######
#  a SmoothCov object containing: 
#    - p by p matrix of sample cov surface estimation on observed grid
#    - NULL for all other entires
##########################################################################

GetCovDense <- function(ymat, mu, optns){
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Sample Covariance is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  n = nrow(ymat)
  K = cov(ymat, use = 'pairwise.complete.obs') # sample variance using non-missing data
  K = 0.5 * (K + t(K))                         # ensure that K is symmetric

  ret = list('rawCov' = NULL, 'smoothCov' = K, 'bwCov' = NULL,
   'sigma2' = NULL, outGrid = NULL)
  class(ret) = "SmoothCov"
  # Garbage Collection
  gc()
  return(ret)
}
