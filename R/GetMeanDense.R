# This function obtains the cross sectional mean function at observed grid
# for dense regular functional data

######
# Input:
######  
#  ymat: matrix of dense regular functional data
#  optns: options for FPCA function
######
# Output: 
######
#  a SMC object containing:
#    - mu: p-dim vector of mean function estimation, i.e. on observed grid
#    - NULL for other entries
##########################################################################

GetMeanDense <- function(ymat, optns){
  # Check optns
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Cross sectional mean is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  mu = colMeans(ymat, na.rm = TRUE) # use non-missing data only
  ret = list('mu' = mu, 'muDense' = NULL, 'mu_bw' = NULL)
  class(ret) = "SMC"
  if(any(is.na(mu))){
      stop('The cross sectional mean is appears to have NaN! Consider setting your dataType to \'Sparse\' manually')
  }
  # Garbage Collection
  gc()
  return(ret)
}
