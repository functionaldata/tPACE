#' Perform FPCA on the functional data 'y' recorderd over 'tt'. Using the options specified in 'p'
#' 
#' @param y is an n-by-1 list of vectors
#' @param tt is an n-by-1 list of vectors
#' @param p is an options structure
#' @return FPCA class model
#' @examples
#' 1 + 3



FPCA = function(y, tt, p = SetOptions()){
  
  # FPCA checks the data validity for the PCA function. 
  if( CheckData(y,tt) ){
    cat('FPCA has stopped.')
    return(TRUE);
  }  
  
  # FPCA checks the options validity for the PCA function. 
  if( CheckOptions(p,num_of_curves) ){
    cat('FPCA has stopped.')
    return(TRUE);
  }

   
  # Conduct FPCA
  #X = PCA(y,tt, bwmu = p$bwmu, bwmu_gcv = p$bwmu_gcv, bwxcov = p$bwxcov,
  #    bwxcov_gcv = p$bwxcov_gcv, ntest1 = p$ntest1, ngrid1 = p$ngrid1, 
  #    selection_k = p$selection_k, FVE_threshold = p$FVE_threshold,
  #    maxk = p$maxk, control = p$control, regular = p$regular, 
  #    error = p$error, ngrid = p$ngrid, method = p$method, 
  #    shrink = p$shrink, newdata = p$newdata, kernel = p$kernel, 
  #    numBins = p$numBins, yname = p$yname, screePlot = p$screePlot, 
  #    designPlot = p$designPlot, rho = p$rho, verbose = p$verbose)
  #
  return(X); 
}
