# sparsify samp
# samp: a matrix of samples, with rows containing the samples
# pts: a vector of grid points, should be from 0 to 1
# sparsity: a vector of integers. The number of observation will be uniform distribution on sparsify.
# aggressive: sparsify in "aggressive" manner making sure that near-by readings are excluded
sparsify <- function(samp, pts, sparsity, aggressive = FALSE) {
    if (length(sparsity) == 1)
            sparsity <- c(sparsity, sparsity) # avoid scaler case
    if(!aggressive){
      indEach <- lapply(1:nrow(samp), function(x) 
        sort(sample(ncol(samp), sample(sparsity, 1))))
    } else {
      indEach <- lapply(1:nrow(samp), function(x)
        remotesampling(ncol(samp), sparsity) )
    }
    tList <- lapply(indEach, function(x) pts[x])
    yList <- lapply(1:length(indEach), function(x) {
        ind <- indEach[[x]]
        y <- samp[x, ind]
        return(y)
    })
   
    return(list(tList=tList, yList=yList))
}



remotesampling <- function(N,s){
  onesamp = sort( sample( N, sample( s, 1))) 
  threshold = (1/length(onesamp))^(1.5) * N 
  while( min(diff(onesamp)) < threshold ){
    onesamp = sort( sample( N, length(onesamp)))
  }
  return(onesamp)
}



