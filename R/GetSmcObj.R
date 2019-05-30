GetSmcObj=function(Ly, Lt, optns = list()){
  CheckData(Ly,Lt)
  
  # Force the data to be list of numeric members and handle NA's
  #Ly <- lapply(Ly, as.numeric) 
  #Lt <- lapply(Lt, as.numeric)
  #Lt <- lapply(Lt, signif, 14)
  #inputData <- list(Ly=Ly, Lt=Lt);
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  optns = SetOptions(Ly, Lt, optns);
  
  # Check the options validity 
  numOfCurves = length(Ly);
  CheckOptions(Lt, optns,numOfCurves)
  
  # Bin the data
  if ( optns$useBinnedData != 'OFF'){ 
    BinnedDataset <- GetBinnedDataset(Ly,Lt,optns)
    Ly = BinnedDataset$newy;
    Lt = BinnedDataset$newt; 
    optns[['nRegGrid']] <- min(optns[['nRegGrid']],
                               BinnedDataset[['numBins']])
    inputData$Ly <- Ly
    inputData$Lt <- Lt
  }
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new
  # data
  # regGrid: the grid of time points for which the smoothed covariance
  # surface assumes values
  # cutRegGrid: truncated grid specified by optns$outPercent for the cov
  # functions
  obsGrid = sort(unique( c(unlist(Lt))));
  regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                          buff & 
                          regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                          buff]
  
  ymat <- List2Mat(Ly, Lt)
  
  ## Mean function
  # If the user provided a mean function use it
  userMu <- optns$userMu
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    smcObj <- GetUserMeanCurve(optns, obsGrid, regGrid, buff)
    smcObj$muDense = ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
  } else if (optns$methodMuCovEst == 'smooth') { # smooth mean
    smcObj = GetSmoothedMeanCurve(Ly, Lt, obsGrid, regGrid, optns)
  } else if (optns$methodMuCovEst == 'cross-sectional') { # cross-sectional mean
    smcObj = GetMeanDense(ymat, obsGrid, optns)
  }
  return(smcObj)
}
