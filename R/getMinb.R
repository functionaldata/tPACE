# In stead of the getMinb.m functionality we can garantee minimal number of neighboring points in here.
getMinb <- function(rcov, out1, regular=rcov$regular, npoly=1, minUniqPts=3, minPts=6) {
    # browser()
    if (regular == 'Sparse') {
        dstar <- minb(out1, 2 + npoly); # rough 1D initial value 
    # get count matrix: TODO rewrite this using rcov$count
        countRes <- getCount(rcov$tpairn)
        distMat <- as.matrix(dist(countRes[, 1:2]))
        # find the bandwidth such that there are at least minUniqPts unique points and minPts points in the 2D window
        bothBW <- sapply(1:ncol(distMat), function(j) {
            # browser()
            x <- distMat[, j]
            ordNeighbors <- order(x)[1:minPts]
            bwNeighbors <- x[ordNeighbors]
            minUniqPtsBW <- bwNeighbors[minPts]
            countNeighbors <- countRes[ordNeighbors, 3]
            minPtsBW <- bwNeighbors[which(cumsum(countNeighbors) >= minPts)[1]]
            return(c(uniqbw=minUniqPtsBW, bw=minPtsBW))
        })
        # browser()
        dstar <- max(dstar, bothBW)
    } else if (regular == 'RegularWithMV') {
        dstar <- minb(out1, 1 + npoly) * 2;
    } else if (regular == 'Dense') {
        dstar = minb(out1, 2 + npoly) * 1.5;
    } 
    
    return(dstar)
}
