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

GetMeanDense <- function(ymat, obsGrid, optns, y = NULL, t = NULL){
  # Check optns
  if(!(optns$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Cross sectional mean is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  if( optns$muCovEstMethod == 'cross-sectional' ){
    if ( is.null(optns$userMu) ){
      mu = colMeans(ymat, na.rm = TRUE) # use non-missing data only
    } else {
      mu = spline(optns$userMu$t, optns$userMu$mu, xout= obsGrid)$y;  
    }
    ret = list('mu' = mu, 'muDense' = NULL, 'mu_bw' = NULL)
  } else if(optns$muCovEstMethod == 'smooth'){
     obsGrid = sort(unique( c(unlist(t))));
     regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
    return( GetSmoothedMeanCurve(y, t, obsGrid, regGrid, optns) )
  } else {
    stop('optns$muCovEstMethod is unknown!\n')
  }
  class(ret) = "SMC"
  if(any(is.na(mu))){
      stop('The cross sectional mean is appears to have NaN! Consider setting your dataType to \'Sparse\' manually')
  }
  # Garbage Collection
  gc()
  return(ret)
}
