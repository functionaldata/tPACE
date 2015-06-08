#' Perform FPCA on the functional data 'y' recorderd over 'tt'. Using the options specified in 'p'
#' 
#' @param y is an n-by-1 list of vectors
#' @param t is an n-by-1 list of vectors
#' @param p is an options structure
#' @return FPCA class model
#' @examples
#' 1 + 3

FPCA = function(y, t, optns = CreateOptions()){
  
  # FPCA checks the data validity for the PCA function. 
  if( CheckData(y,t) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }  
 

  # FPCA sets the options structure that are still NULL
  optns = SetOptions(y, t, optns);

  
  # FPCA checks the options validity for the PCA function. 
  numOfCurves = length(y);
  if( CheckOptions(t, optns,numOfCurves) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }


  # Bin the data (potentially):
  if ( optns$useBinnedData != 'OFF'){ 
      BinnedDataset <- GetBinnedDataset(y,t,optns)
      y = BinnedDataset$newy;
      t = BinnedDataset$newt; 
  }

  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new data
  # regGrid: the grid of time points for which the smoothed covariance surface assumes values
  obsGrid = sort(unique( c(unlist(t), optns$newdata)));
  regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);


  # Get the smoothed mean curve
  smcObj = GetSmoothedMeanCurve(y, t, obsGrid, regGrid, optns)


  # Get the smoothed covariance surface
  # mu: the smoothed mean curve evaluated at times 'obsGrid'
  mu = smcObj$mu
  scsObj = GetSmoothedCovarSurface(y, t, mu, obsGrid, regGrid, optns) 

  # internal working regular grid. This may be a truncated version of the regGrid.
  workGrid <- scsObj$outGrid

  # Get the results for the eigen-analysis
  eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, workGrid, optns)

  # truncated obsGrid, and observations. Empty observation due to truncation has length 0.
  if (!all.equal(optns$outPercent, c(0, 1))) {
    buff <- .Machine$double.eps * 10
    obsGrid <- obsGrid[obsGrid >= min(workGrid) - buff &
                            obsGrid <= max(workGrid) + buff]
    tmp <- TruncateObs(y, t, obsGrid)
    y <- tmp$y
    t <- tmp$y
  }

  # convert things
  # convert phi and fittedCov to obsGrid.
  phiObs <- ConvertSupport(workGrid, obsGrid, eigObj$phi)
  CovObs <- ConvertSupport(workGrid, obsGrid, eigObj$fittedCov)

  # Get scores

  # Make the return object X
  X <- list( 'sigma' = scsObj$sigma, 'eigVal' = eigObj$eigVal, 'eigFunc' = eigObj$eigFunc, 'mu' = smcObj$mu, 
             'smoothedCov' = scsObj$smoothCov, 'fittedCov' = eigObj$fittedCov)

  return(X); 
}
