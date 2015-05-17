gcvlwls2d <- function(tt, ngrid=NULL, regular=rcov$regular, error=rcov$error, kern, rcov, h0=getMinb(rcov, regular=rcov$regular, out1=sort(unique(unlist(tt)))), verbose=TRUE) {
# TODO: Consider rewrite the function into a more general form. The current implementation is too specific and works only for smoothing raw covariance.

# Returns: a list of length 2, containing the optimal bandwidth and the gcv score.
# tt: input time points. 
# ngrid: I think this should not be used in the gcv function.

# This function computes the optimal bandwidth choice for the covariance surface. 
# function use GCV method by pooling the longitudinal data together. 
# verbose is unused for now
# this is incompatible with PACE because the GCV is calculated in a different way.
  
t <- unlist(tt)
r <- diff(range(t)) * sqrt(2) # sqrt(2) because the window is circular.

if (kern == 'gauss') {
    if (is.null(h0))
        stop('Not implemented')
    
    h0 = h0 * 0.2;
}

if (is.null(h0))
    stop('the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n')

# TODO: improve this initial choice. The initial h0 may not be that good. The search scheme may be too rough if there are a lot number of observations.
h0 <- min(h0, r/4)
if (h0 < r/4) {    
    q <- (r / (4 * h0)) ^ (1/9)
} else if (h0 < r/2) {
    q <- (r / (2 * h0)) ^ (1/9)
} else if (h0 < r) {
    q <- (r / h0) ^ (1/9)
} else {
    stop('Data is too sparse. The minimal bandwidth is the range of data')
}

bw <- (q ^ (0:9)) * h0 # from h0 to r / 4

opth <- h0

leave <- FALSE
iter <- 0
maxIter <- 1
while (!leave && iter < maxIter) {
    gcvScores <- sapply(bw, getGCVscores, kern=kern, xin=rcov$tpairn, yin=rcov$cxxn, win=as.numeric(rcov$win)) 
    gcvFrame <- data.frame(bw=bw, gcv=gcvScores)
    optInd <- which.min(gcvScores)
    opth <- bw[optInd]
    optgcv <- gcvScores[optInd]
    
    if (opth >= r - 1e-12) {
        leave <- TRUE
        stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
    }
    # else if (opth < r) {
    
    if (optInd != length(bw) && !is.infinite(optgcv))
        leave <- TRUE            
    else if (is.infinite(optgcv)) {
        warnings('Data is too sparse, retry with larger bandwidths!')
        h0 <- bw[10] * 1.01
    } else if (opth == bw[length(bw)]) {
        warning('Optimal bandwidth not found in the candidate bandwidths. Retry with larger bandwidths')
        h0 <- bw[9] 
    }
    
    if (!leave) {
        newr <- seq(0.5, 1, by=0.05) * r # ??? this can be quite slow
        ind <- which(newr > h0)[1]
        q <- (newr[ind] / h0) ^ (1/9)
        bw <- q ^ (0:9) * h0
        if (verbose) {
            cat('New bwxcov candidates:\n')
            print(bw)
        }

        iter <- iter + 1
    }
    # }
        
}

return(list(h=opth, gcv=optgcv))

}


getGCVscores <- function(bw, ...) {
# ...: passed on to lwls2d
    
    fit <- lwls2d(bw, ..., returnFit=TRUE)
    return(gcv(fit)['gcv'])

}


