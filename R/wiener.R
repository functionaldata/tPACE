# TODO: Roxygen documentation
#' @export

# A test on standard Wiener process (brownian motion)
# n: sample sizeDiss
# pts: a vector of grid points, should be from 0 to 1
# K: number of components
# sparsify: a vector of integers. The number of observation will be uniform distribution on sparsify.
wiener <- function(n=1, pts=seq(0, 1, length=50), sparsify=NULL, K=50) {
    # Simulate n standard Wiener process on [0, 1] observed at pts. K is the number of components used in the simulation.
    # Each row is a sample.
    pts <- as.matrix(pts)
    if (dim(pts)[1] < dim(pts)[2])
        pts <- t(pts)
        
    basis <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
    samp <- t(basis %*% diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
    
    if (!is.null(sparsify)) {
        samp <- sparsify(samp, pts, sparsify)
    }
    
    return(samp)

}

## sparsify samp
## samp: a matrix of samples, with rows containing the samples
## pts: a vector of grid points, should be from 0 to 1
## sparsity: a vector of integers. The number of observation will be uniform distribution on sparsify.
#sparsify <- function(samp, pts, sparsity) {
#    if (length(sparsity) == 1)
#            sparsity <- c(sparsity, sparsity) # avoid scaler case
#    
#    indEach <- lapply(1:nrow(samp), function(x) 
#        sort(sample(ncol(samp), sample(sparsity, 1))))
#    tList <- lapply(indEach, function(x) pts[x])
#    yList <- lapply(1:length(indEach), function(x) {
#        ind <- indEach[[x]]
#        y <- samp[x, ind]
#        return(y)
#    })
#   
#    return(list(tList=tList, yList=yList))
#}

