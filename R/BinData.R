BinData = function(y,t,regular,verbose,numBins){
  
  # Bin the data 'y'
  # y : n-by-1 list of vectors
  # t : n-by-1 list of vectors 
  # regular : indicator about structure of the data 
  #   (dense (2), or  regular data with missing values (1) or sparse (0))
  # verbose gives warning messages
  # numBins: number of bins (if set)

  if( !(  any(regular == c(1,2,0)) )){  
    cat("Error: BinData is aborted because the argument: regular is invalid!\n");  
    return(TRUE);   
  }

  BinDataOutput <- list( newy <- NULL, newt <- NULL);
  
  n = length(t);
  ni = sapply(FUN= length,t);
  
  if (is.null(regular)){
    regular = 0;
  }
  
  if (is.null(verbose)){
    verbose = TRUE;
  }
  
  if (regular == 0){
    m = median(ni)
  } else {
    m = max(ni);
  }
  
  if (is.null(numBins)){
    numBins = GetBinNum(n,m,regular,verbose)
  }else if( numBins <1){
    cat("Warning: number of bins must be positive integer! Reset to default number of bins now!\n")
    numBins = GetBinNum(n,m,regular,verbose)
  }

  if (is.null(numBins)){
    return( BinDataOutput )
  }
     
  numBins = ceiling(numBins);
  
  tt = unlist(t);  
  a0 = min(tt);
  b0 = max(tt);
  
  if (verbose){
    cat("We start binning the data!\n")
  } 
  
  for (i in 1:n){
    res = GetBinnedCurve(t[[i]], y[[i]], numBins, TRUE, TRUE, a0, b0);
    BinDataOutput$newt[i] = res$midpoint;   
    BinDataOutput$newy[i] = res$newy;      
  }
  if (verbose){ 
    cat("We finished binning the data in the time domain [", a0,",", b0, "] using", numBins, "bins.\n" )
  } 
  
  return( BinDataOutput )
}
