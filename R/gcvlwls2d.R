gcvlwls2d <- function(tt, ngrid=NULL, regular=rcov$regular, error=rcov$error, kern, rcov, h0=NULL, useBins=FALSE, verbose=TRUE, CV=FALSE) {
# TODO: Consider rewrite the function into a more general form. The current implementation is too specific and works only for smoothing raw covariance.
# TODO: implement useBins 

# Returns: a list of length 2, containing the optimal bandwidth and the gcv score.
# tt: input time points. 
# ngrid: I think this should not be used in the gcv function.
# CV: whether to use CV rather than GCV. Supported values are '[k]fold', where '[k]' is a integer.

# This function computes the optimal bandwidth choice for the covariance surface. 
# function use GCV method by pooling the longitudinal data together. 
# verbose is unused for now
# this is incompatible with PACE because the GCV is calculated in a different way.
  
  t <- unlist(tt)
  r <- diff(range(t)) * sqrt(2) # sqrt(2) because the window is circular.

  minBW <- getMinb(rcov, regular=rcov$regular, out1=sort(unique(unlist(tt))))

  if (missing(h0)) {
    h0 <- minBW
  }
  
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
  maxIter <- 1 # ??
  
  if (CV != FALSE) {
    useKfold <- regexpr('fold', CV)
    if (useKfold != -1) {
      fold <- as.integer(substr(CV, 1, useKfold - 1))
      partition <- caret::createFolds(rcov$cxxn, k=fold)
    }
  }
  
  while (!leave && iter < maxIter) {
    if (CV == FALSE) {
      Scores <- sapply(bw, getGCVscores, kern=kern, xin=rcov$tpairn, yin=rcov$cxxn, win=as.numeric(rcov$win)) 
    } else {
      Scores <- sapply(bw, getCVscores, partition=partition, kern=kern, xin=rcov$tpairn, yin=rcov$cxxn, win=as.numeric(rcov$win)) 
    }
    optInd <- which.min(Scores)
    opth <- bw[optInd]
    optgcv <- Scores[optInd]
    
    if (opth >= r - 1e-12) {
      leave <- TRUE
      stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
    }
    # else if (opth < r) {
      
      if (optInd != length(bw) && !is.infinite(optgcv))
      leave <- TRUE            
      else if (is.infinite(optgcv)) {
        if (verbose)
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

  ret <- list(h=opth, gcv=optgcv, minBW=minBW)
  if (CV != FALSE)
    names(ret)[2] <- 'cv'
  
  return(ret)

}


getGCVscores <- function(bw, ...) {
# ...: passed on to lwls2d
  
  fit <- lwls2d(bw, ..., returnFit=TRUE)
  return(gcv(fit, maxk=fit$frame$maxk)['gcv'])

}


# k-fold CV
# partition: a list of testset observation indices, returned by caret::createFolds
# ...: passed on to lwls2d
getCVscores <- function(bw, partition, xin, yin, ...) {
  n <- length(yin)
  
  # browser()
  cvSubSum <- sapply(partition, function(testSet) {
    # browser()
    fit <- lwls2d(bw, xin=xin, yin=yin, subset=-testSet, ..., returnFit=TRUE)
    cvSubSum <- (yin[testSet] - predict(fit, newdata=xin[testSet, ]))^2
    return(sum(cvSubSum))
  })
  
  return(sum(cvSubSum))
}
