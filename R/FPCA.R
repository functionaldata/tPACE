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


  if(optns$dataType == 'Dense' || optns$dataType == 'DenseWithMV'){
    # Only applicable to Dense and Regular functional data,
    # or Dense data with Missing Values 
    
    # cross sectional mean and sample covariance for dense case
    # assume no measurement error.
    ymat = List2Mat(y)

    # Define time grids
    obsGrid = sort(unique(unlist(t)))
    regGrid = obsGrid
    workGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid)

    # get cross sectional mean and sample cov
    smcObj = GetMeanDense(ymat, optns)
    mu = smcObj$mu
    smcObj$muDense = ConvertSupport(obsGrid, workGrid, mu = mu)
    scsObj = GetCovDense(ymat, mu, optns)
    scsObj$smoothCov = ConvertSupport(obsGrid, workGrid, Cov = scsObj$smoothCov)

  } else if(optns$dataType == 'Sparse'){
    # For Sparse case
    
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
    scsObj = GetSmoothedCovarSurface(y, t, mu, obsGrid, regGrid, optns, optns$useBins) 
    sigma2 <- scsObj$sigma2

    # workGrid: possibly truncated version of the regGrid; truncation would occur during smoothing
    workGrid <- scsObj$outGrid
  } else {
    stop('not implemented for dataType = "p>>n" yet!')
  }

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

  # convert phi and fittedCov to obsGrid.
  phiObs <- ConvertSupport(workGrid, obsGrid, phi=eigObj$phi)
  CovObs <- ConvertSupport(workGrid, obsGrid, Cov=eigObj$fittedCov)

  # Get scores  
  if (optns$method == 'CE') {
    if (optns$rho != 'no') {
      rho <- GetRho(y, t, optns, mu, obsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      sigma2 <- rho
    }

    scoresObj <- GetCEScores(y, t, optns, mu, obsGrid, CovObs, eigObj$lambda, phiObs, sigma2)

  } else if (optns$method == 'IN') {
    scoresObj <- GetINScores(ymat, t, optns, mu, eigObj$lambda, phiObs)
  }

  # Make the return object by MakeResultFPCA
  ret <- MakeResultFPCA(optns, smcObj, mu, scsObj, eigObj,
  scoresObj, obsGrid, workGrid)

  return(ret); 
}
