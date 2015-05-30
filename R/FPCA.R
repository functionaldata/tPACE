#' Perform FPCA on the functional data 'y' recorderd over 'tt'. Using the options specified in 'p'
#' 
#' @param y is an n-by-1 list of vectors
#' @param t is an n-by-1 list of vectors
#' @param p is an options structure
#' @return FPCA class model
#' @examples
#' 1 + 3



FPCA = function(y, t, p = NULL){
  
  # FPCA checks the data validity for the PCA function. 
  if( CheckData(y,t) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }  
 
  # FPCA sets the options structure that are still NULL
  p = SetOptions(y, t, p);
  
  # FPCA checks the options validity for the PCA function. 
  numOfCurves = length(y);
  if( CheckOptions(t, p,numOfCurves) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }


  # Bin the data (potentially):
  if ( p$use_binned_data != 'OFF'){ 
      BinnedDataset <- GetBinnedDataset(y,t,p)
      y = BinnedDataset$newy;
      t = BinnedDataset$newt; 
  }

  # Generate basic grids:
  # out1:  the unique sorted pooled time points of the sample and the new data
  # out21: the grid of time points for which the smoothed covariance surface assumes values
  out1 = sort(unique( c(unlist(t), p$newdata)));
  out21 = seq(min(out1), max(out1),length.out = p$ngrid);


  # Get the smoothed mean curve
  smcObj = GetSmoothedMeanCurve(y, t, out1, out21, p)


  # Get the smoothed covariance surface
  # mu: the smoothed mean curve evaluated at times 'out1'
  mu = smcObj$mu
  scsObj = GetSmoothedCovarSurface(y, t, out1, mu, p) 


  # Get the results for the eigen-analysis
  eigObj = GetEigenAnalysisResults(y,t, scsObj, out21, p)


  # Make the return object X
  X <- list( 'sigma' = scsObj$sigma, 'eigVal' = eigObj$eigVal, 'eigFunc' = eigObj$eigFunc, 'mu' = smcObj$mu, 
             'smoothedCov' = scsObj$xcov, 'fittedCov' = scsObj$xcovfit)

  return(X); 
}
