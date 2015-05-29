GetBinnedDataset <- function (y, t, p){
   
  # Bin the data 'y'
  # y : n-by-1 list of vectors
  # t : n-by-1 list of vectors 

  BinDataOutput <- list( newy <- NULL, newt <- NULL);

  regular = p$regular;
  verbose = p$verbose;  
  numBins = p$numBins;
  tt = unlist(t);  
  a0 = min(tt);
  b0 = max(tt);

  n = length(t);
  ni = sapply(FUN= length,t);
  
  if (regular == 'Sparse'){
    m = median(ni)
  } else {
    m = max(ni);
  }
  
  # Determine the number of bins automatically if numBins is null
  if (is.null(numBins) && p$use_binned_data =='AUTO'){
    numBins = GetBinNum(n,m,regular,verbose)
    # and if it is still NULL return the unbinned data
    if (is.null(numBins)){
      BinDataOutput$newt = t;
      BinDataOutput$newy = y;
      return( BinDataOutput )
    }
  }
  # otherwise use the one provided by the user (ceiled)
  numBins = ceiling(numBins);

  for (i in 1:n){
    res = GetBinnedCurve(t[[i]], y[[i]], numBins, TRUE, TRUE, a0, b0);
    BinDataOutput$newt[i] = res$midpoint;   
    BinDataOutput$newy[i] = res$newy;      
  }
     
  result <- list( 't' = BinDataOutput$newt, 'y' = BinDataOutput$newy);
  return(result)
}


