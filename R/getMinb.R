# In stead of the getMinb.m functionality we can garantee minimal number of neighboring points in here.
# TODO: distMat is memory inefficient.
getMinb <- function(rcov, obsGrid, dataType=rcov$dataType, npoly=1, minUniqPts=3, minPts=6) {
  # browser()
  if (dataType == 'Sparse') {
    dstar <- minb(obsGrid, 2 + npoly) # rough 1D initial value 
    # get count matrix: TODO rewrite this using rcov$count
    if (class(rcov) == 'RawCov') {
      countRes <- getCount(rcov$tpairn)
      count <- countRes[, 3]
      distMat <- as.matrix(dist(countRes[, 1:2]))
    } else if (class(rcov) == 'BinnedRawCov') {
      count <- rcov$count
      distMat <- as.matrix(dist(rcov$tPairs))
    }
    # find the bandwidth such that there are at least minUniqPts unique points and minPts points in the 2D window
    bothBW <- sapply(1:ncol(distMat), function(j) {
      # browser()
      x <- distMat[, j]
      ordNeighbors <- order(x)[1:minPts]
      bwNeighbors <- x[ordNeighbors]
      minUniqPtsBW <- bwNeighbors[minPts]
      countNeighbors <- count[ordNeighbors]
      minPtsBW <- bwNeighbors[which(cumsum(countNeighbors) >= minPts)[1]]
      return(c(uniqbw=minUniqPtsBW, bw=minPtsBW))
    })
    # browser()
    dstar <- max(dstar, bothBW)
  } else if (dataType == 'RegularWithMV') {
    dstar <- minb(obsGrid, 1 + npoly) * 2;
  } else if (dataType == 'Dense') {
    dstar = minb(obsGrid, 2 + npoly) * 1.5;
  } 
  
  return(dstar)
}
