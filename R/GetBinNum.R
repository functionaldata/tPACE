GetBinNum = function(n, m, regular, verbose ){
  
  # Get the number of bins 
  # n : number of curves
  # m : median or max value of number of time-points
  # regular : indicator about structure of the data 
  #     (dense (2), or  regular data with missing values (1) or sparse (0))
  # verbose : outpit diagnostics/progress
  
  numBin = NULL;
  if (m <= 20){
    if (regular ==0){
      str = 'Median of ni';
    } else {
      str = 'Maximum of ni';
    }
    if (verbose){
      cat(str, 'is no more than 20! No binning is performed!\n');    
    }
    return(NULL)
  }

  
  if (m >400){
    numBin = 400;
  }
  
  if (n > 5000){
    mstar = max(20,(((5000-n)*19)/2250)+400);
    if (mstar < m){
      numBin = ceiling(mstar);
    } else {
      if (verbose){
        cat('No binning is needed!\n');    
      }
      return(NULL)
    }
  }
  
  if( verbose && is.null(numBins) ) {   
    cat('No binning is needed!\n');
  }

  return(numBin)
}
