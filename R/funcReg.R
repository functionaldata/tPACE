## Concurrent Functional regression by 2D smoothing method.
# vars: a list of input functional/scaler covariates. Each field corresponds to a functional (a list) or scaler (a vector) covariate. The last entry is assumed to be the response if no entry is names 'Y'. If a field corresponds to a functional covariate, it should have two fields: 'tList', a list of time points, and 'yList', a list of function values.
# bw: bandwidth used.
# Tout: output time points.
# kern: kernel used.
# measurementError: whether measurement errors on the functional observations should be assumed. If TRUE the diagonal raw covariance will be removed when smoothing.
# diag1D: A string specifying whether to use 1D smoothing for the diagonal line of the covariance. 'none': don't use 1D smoothing; 'cross': use 1D only for cross-covariances; 'all': use 1D for both auto- and cross-covariances.
# returnCov: whether to return the covariance surfaces, which is a four dimensional array. The first two dimensions correspond to Tout, and the last two correspond to the covariates and the response, i.e. (i, j, k, l) entry being Cov(X_k(t_i), X_l(t_j))
mvConReg <- function(vars, bw, Tout, kern='gauss', measurementError=TRUE, diag1D='none', returnCov=TRUE) {
  
  n <- lengthVars(vars)
  p <- length(vars) - 1
  if (p == 0)
    stop('Too few covariates.')
  
  if (is.null(names(vars)))
    names(vars) <- c(paste0('X', seq_len(length(vars) - 1)), 'Y')
  
  if ('Y' %in% names(vars)) {
    vars <- c(vars[names(vars) != 'Y'], vars['Y'])
  } else if (names(vars)[length(vars)] == '') {
    names(vars)[length(vars)] <- 'Y'
  }
  
  Yname <- names(vars)[length(vars)]
  
  # De-mean.
  demeanedRes <- demean(vars, bw, kern)
  vars <- demeanedRes[['xList']]
  muList <- demeanedRes[['muList']]
  
  allCov <- MvCov(vars, bw, Tout, kern, measurementError, center=FALSE, diag1D)
  beta <- sapply(seq_len(dim(allCov)[1]), function(i) {
    tmpCov <- allCov[i, i, , ]
    beta_ti <- qr.solve(tmpCov[1:p, 1:p], tmpCov[1:p, p + 1])
    beta_ti
  })
  if (is.null(nrow(beta)))
    beta <- matrix(beta, 1)
  rownames(beta) <- names(vars)[-length(vars)]
  
  # coefficient of determination: 
  #   R2 = cov(X, Y)' var(X)^{-1} cov(X, Y) / var(Y)
  R2 <- sapply(seq_len(dim(allCov)[1]), function(i) {
    tmpCov <- allCov[i, i, , ]
    tmpCov[p + 1, 1:p, drop=FALSE] %*% beta[, i, drop=FALSE] / tmpCov[p + 1, p + 1]
  })
  
  muBeta <- sapply(seq_len(p), function(j) {
    if (!is.function(muList[[j]])) { # scaler mean
      beta[j, ] * rep(muList[[j]], length(Tout))
    } else { # functional mean
      beta[j, ] * muList[[j]](Tout)
    }
  })
  alpha <- muList[[Yname]](Tout) - colSums(t(muBeta))
  
  res <- list(beta=beta, alpha=alpha, Tout=Tout, cov=allCov, R2=R2, n=n)
  if (!returnCov)
    res[['cov']] <- NULL
  res
}

demean <- function(vars, bw, kern) {
  tmp <- lapply(vars, function(x) {
    if (is.numeric(x)) { # scaler
      xmu <- mean(x)
      x <- x - xmu
    } else if (is.list(x)) { # functional
      Tin <- sort(unique(unlist(x[['tList']])))
      xmu <- GetSmoothedMeanCurve(x[['yList']], x[['tList']], Tin, Tin[1],
                                  list(userBwMu=bw, kernel=kern))[['mu']]
      muFun <- approxfun(Tin, xmu)
      x[['yList']] <- lapply(1:length(x[['yList']]), function(i)
        x[['yList']][[i]]- muFun(x[['tList']][[i]]))
      xmu <- muFun
    }
    
    list(x=x, mu=xmu)
  })
  
  xList <- lapply(tmp, `[[`, 'x')
  muList <- lapply(tmp, `[[`, 'mu')
  
  list(xList = xList, muList = muList)
}

## Multivariate function/scaler covariance.
# INPUTS: same as mvConReg
# Output: a 4-D array containing the covariances. The first two dimensions corresponds to 
# time s and t, and the last two dimensions correspond to the variables taken covariance upon.
MvCov <- function(vars, bw, Tout, kern, measurementError=TRUE, center=TRUE, diag1D='none') {
  if (!is.list(vars) || length(vars) < 1)
    stop('`vars` needs to be a list of length >= 1')
  
  isFuncVars <- sapply(vars, is.list)
  p <- length(isFuncVars)
  pFunc <- sum(isFuncVars)
  pScaler <- sum(!isFuncVars)
  
  if (any(isFuncVars)) {
    tAll <- do.call(c, lapply(vars[isFuncVars], function(x) unlist(x[['tList']])))
    Tin <- sort(unique(tAll))
    
    if (missing(Tout))
      Tout <- Tin
    lenTout <- length(Tout)
    
  } else {
    stop('No functional observation found')
  }
  
  # First two dimensions are for s, t, and the last two dimensions are for matrix of # random variables.
  res <- array(NA, c(lenTout, lenTout, p, p))
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      if (j <= i) {
        use1D <- diag1D == 'all' || ( diag1D == 'cross' && j != i )
        covRes <- uniCov(vars[[i]], vars[[j]], bw, Tout, kern, 
                         rmDiag = (i == j) && measurementError, 
                         center, use1D)
        if (attr(covRes, 'covType') %in% c('FF', 'SS'))
          res[, , i, j] <- covRes
        else {
          if (nrow(covRes) == 1)   # cov(scaler, function)
            res[, , i, j] <- matrix(covRes, lenTout, lenTout, byrow=TRUE)
          else                     # cov(function, scaler)
            res[, , i, j] <- matrix(covRes, lenTout, lenTout, byrow=FALSE)
        }
      } else { # fill up the symmetric cov(y, x)
        res[, , i, j] <- t(res[, , j, i])
      }
    }
  }
  
  return(res)
}

## Univariate function/scaler covariance.
# rmDiag: whether to remove the diagonal of the raw covariance. Ignored if 1D smoother is used.
# center: whether to center the covariates before calcuate covariance.
# use1D: whether to use 1D smoothing for estimating the diagonal covariance.
uniCov <- function(X, Y, bw, Tout, kern='gauss', rmDiag=FALSE, center=TRUE, use1D=FALSE) {
  flagScalerFunc <- FALSE
  # Force X to be a function in the scaler-function case.
  if (!is.list(X) && is.list(Y)) {
    flagScalerFunc <- TRUE
    tmp <- X
    X <- Y
    Y <- tmp
  }
  
  # Scaler-scaler
  if (!is.list(X) && !is.list(Y)) {
    res <- cov(X, Y)
    attr(res, 'covType') <- 'SS'
    
    # Scaler-function    
  } else if (is.list(X) && !is.list(Y)) {
    Tin <- sort(unique(unlist(X[['tList']])))
    if (center) {
      Xmu <- GetSmoothedMeanCurve(X[['yList']], X[['tList']], Tin, Tin[1], list(userBwMu=bw, kernel=kern))[['mu']]
      Ymu <- mean(Y)
    } else {
      Xmu <- rep(0, length(Tin))
      Ymu <- 0
    }
    res <- GetCrCovYZ(bw, Y, Ymu, X[['yList']], X[['tList']], Xmu, Tin, kern)[['smoothedCC']]
    res <- as.matrix(ConvertSupport(Tin, Tout, mu=res))
    if (flagScalerFunc) 
      res <- t(res)
    
    attr(res, 'covType') <- 'FS'
    
    # function-function  
  } else {
    TinX <- sort(unique(unlist(X[['tList']])))
    TinY <- sort(unique(unlist(Y[['tList']])))
    nTout <- length(Tout)
    if (center) {
      if (min(TinX) > min(Tout) || min(TinY) > min(Tout) || 
          max(TinY) < max(Tout) || max(TinX) < max(Tout))
        stop('Observation time points coverage too low')
      
      Xmu <- GetSmoothedMeanCurve(X[['yList']], X[['tList']], TinX, TinX[1],
                                  list(userBwMu=bw, kernel=kern))[['mu']]
      Ymu <- GetSmoothedMeanCurve(Y[['yList']], Y[['tList']], TinY, TinY[1],
                                  list(userBwMu=bw, kernel=kern))[['mu']]
    } else {
      Xmu <- rep(0, length(TinX))
      Ymu <- rep(0, length(TinY))
    }
    names(Xmu) <- TinX
    names(Ymu) <- TinY
    
    if (use1D) {
      Xvec <- unlist(X[['yList']])
      Yvec <- unlist(Y[['yList']])
      tvecX <- unlist(X[['tList']])
      tvecY <- unlist(Y[['tList']])
      if (!identical(tvecX, tvecY)){
        stop('Cannot use 1D covariance smoothing if the observation time points for X and Y are different')
      }
      
      ord <- order(tvecX)
      tvecX <- tvecX[ord]
      Xvec <- Xvec[ord]
      Yvec <- Yvec[ord]
      Xcent <- Xvec - Xmu[as.character(tvecX)]
      Ycent <- Yvec - Ymu[as.character(tvecX)]
      covXY <- Lwls1D(bw, kern, npoly=1L, nder=0L, 
                      xin=tvecX, yin=Xcent * Ycent, 
                      win=rep(1, length(tvecX)), xout=Tout)
      res <- matrix(NA, nTout, nTout)
      diag(res) <- covXY
    } else { # use 2D smoothing
      tmp <- GetCrCovYX(bw, bw, X[['yList']], X[['tList']], Xmu, 
                        Y[['yList']], Y[['tList']], Ymu, rmDiag=rmDiag, kern=kern)
      gd <- tmp[['smoothGrid']]
      res <- matrix(interp2lin(gd[, 1], gd[, 2], tmp[['smoothedCC']], rep(Tout, times=nTout), rep(Tout, each=nTout)), nTout, nTout)
    }
    attr(res, 'covType') <- 'FF'
  }
  
  return(res)
}

## Concurrent functional regression by imputation. This does not provide consistent estimates.
## FPCAlist: a list of functional covariates and response. Each field corresponds to a covariate. 
#            The last entry is assumed to be the response if no entry is names 'Y'.
imputeConReg <- function(FPCAlist, Z, Tout) {
  
  if (is.null(names(FPCAlist)))
    names(FPCAlist) <- c(paste0('X', seq_len(length(FPCAlist) - 1)), 'Y')
  
  if ('Y' %in% names(FPCAlist)) {
    Yname <- 'Y'
    FPCAlist <- c(FPCAlist[names(FPCAlist) != 'Y'], FPCAlist['Y'])
  } else 
    Yname <- names(FPCAlist)[length(FPCAlist)]
  
  imputeCurves <- sapply(FPCAlist, function(x) 
    apply(fitted(x), 1, function(fit) 
      approx(x[['workGrid']], fit, Tout)[['y']]),
    simplify='array')
  alphaBeta <- apply(imputeCurves, 1, function(XYt) {
    Yt <- XYt[, ncol(XYt)]
    designMat <- cbind(1, XYt[, -ncol(XYt), drop=FALSE], Z)
    beta_t <- qr.solve(designMat, Yt)
    return(beta_t)
  })
  beta0 <- alphaBeta[1, ]
  beta <- alphaBeta[-1, , drop=FALSE]
  
  return(list(beta0 = beta0, beta = beta, Tout = Tout))
}

## regObj: an object returned by mvConReg.
## vars: a list of input functional/scaler covariates. Each field  can correspond to a covariate. 
#        The last entry is assumed to be the response if no entry is names 'Y'.
summaryConReg <- function(regObj, vars) {
  
}

## subset a list of covariates and responses.
subsetVars <- function(vars, subset) {
  sapply(vars, function(x) {
    if (is.list(x)) {
      sapply(x, `[`, subset, drop=FALSE, simplify=FALSE)
    } else if (is.numeric(x)) {
      x[subset, drop=FALSE]
    } else {
      stop('Cannot subset variable')
    }
  }, simplify=FALSE)
}

## get the number of subjects for a list of covariates and responses.
lengthVars <- function(vars, subset) {
  lenEach <- sapply(vars, function(x) {
    if (is.list(x)) {
      sapply(x, length)
    } else if (is.numeric(x)) {
      length(x)
    } else {
      stop('Cannot subset variable')
    }
  }, simplify=FALSE)
  len <- unique(unlist(lenEach))
  if (length(len) != 1) {
    stop('Length of variables are not the same!')
  }
  
  return(len)
}
