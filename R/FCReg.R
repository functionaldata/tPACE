#' Functional Concurrent Regression by 2D smoothing method.
#' 
#' Functional concurrent regression with dense or sparse functional data for scalar or functional dependent variable. 
#' 
#' @param vars A list of input functional/scaler covariates. Each field corresponds to a functional (a list) or scaler (a vector) covariate. The last entry is assumed to be the response if no entry is names 'Y'. If a field corresponds to a functional covariate, it should have two fields: 'Lt', a list of time points, and 'Ly', a list of function values.
#' @param userBwMu A scalar with bandwidth used for smoothing the mean
#' @param userBwCov A scalar with bandwidth used for smoothing the auto- and cross-covariances
#' @param outGrid A vector with the output time points
#' @param kern Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" (default: "gauss")
#' @param measurementError Indicator measurement errors on the functional observations should be assumed. If TRUE the diagonal raw covariance will be removed when smoothing. (default: TRUE)
#' @param diag1D  A string specifying whether to use 1D smoothing for the diagonal line of the covariance. 
#' 'none': don't use 1D smoothing; 'cross': use 1D only for cross-covariances; 'all': use 1D for both auto- and cross-covariances. (default : 'none')
#' @param useGAM Indicator to use gam smoothing instead of local-linear smoothing (semi-parametric option) (default: FALSE)
#' @param returnCov Indicator to return the covariance surfaces, which is a four dimensional array. The first two dimensions correspond to outGrid
#'  and the last two correspond to the covariates and the response, i.e. (i, j, k, l) entry being Cov(X_k(t_i), X_l(t_j)) (default: FALSE)
#' @param ...  Additional arguments 
#' 
#' @details If measurement error is assumed, the diagonal elements of the raw covariance will be removed. This could result in highly unstable estimate if the design is very sparse, or strong seasonality presents. 
#' @references
#' \cite{Yao, F., Mueller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#'
#' \cite{Senturk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)} 
#' @export
#' @examples 
#' # Y(t) = \beta_0(t) + \beta_1(t) X_1(t) + \beta_2(t) Z_2 + \epsilon
#' 
#' # Settings
#' set.seed(1)
#' n <- 100
#' nGridIn <- 200
#' sparsity <- 5:10 # Sparse data sparsity
#' T <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
#' bw <- 0.1
#' outGrid <- round(seq(min(T), 1, by=0.05), 2)
#' 
#' # Simulate functional data 
#' mu <- T * 2 # mean function for X_1
#' sigma <- 1
#' 
#' beta_0 <- 0
#' beta_1 <- 1
#' beta_2 <- 1
#' 
#' Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
#' X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
#' epsilon <- rnorm(n, sd=sigma)
#' Y <- matrix(NA, n, nGridIn)
#' for (i in seq_len(n)) {
#'   Y[i, ] <- beta_0 + beta_1 * X_1[i, ] + beta_2 * Z[i, 2] + epsilon[i]
#' }
#' 
#' # Sparsify functional data
#' set.seed(1)
#' X_1sp <- Sparsify(X_1, T, sparsity)
#' set.seed(1)
#' Ysp <- Sparsify(Y, T, sparsity)
#' vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
#' withError2D <- FCReg(vars, bw, bw, outGrid)


FCReg <- function(vars, userBwMu, userBwCov, outGrid, kern='gauss', measurementError=TRUE, diag1D='none', useGAM = FALSE, returnCov=TRUE) {
  
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

  # Handle NaN, int to double
  vars[sapply(vars, is.list)] <- lapply(
    vars[sapply(vars, is.list)], 
    function(v) HandleNumericsAndNAN(v[['Ly']], v[['Lt']])
  )
  # outGrid <- as.numeric(outGrid)
  
  # De-mean.
  demeanedRes <- demean(vars, userBwMu, kern)
  vars <- demeanedRes[['xList']]
  muList <- demeanedRes[['muList']]
  
  allCov <- MvCov(vars, userBwCov, outGrid, kern, measurementError, center=FALSE, diag1D)
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
      beta[j, ] * rep(muList[[j]], length(outGrid))
    } else { # functional mean
      beta[j, ] * muList[[j]](outGrid)
    }
  })
  beta0 <- muList[[Yname]](outGrid) - colSums(t(muBeta))
  
  res <- list(beta=beta, beta0 = beta0, outGrid=outGrid, cov=allCov, R2=R2, n=n)
  if (!returnCov)
    res[['cov']] <- NULL
  res
}

demean <- function(vars, userBwMu, kern) {
  tmp <- lapply(vars, function(x) {
    if (is.numeric(x)) { # scaler
      xmu <- mean(x)
      x <- x - xmu
    } else if (is.list(x)) { # functional
      Tin <- sort(unique(unlist(x[['Lt']])))
      xmu <- GetSmoothedMeanCurve(x[['Ly']], x[['Lt']], Tin, Tin[1],
                                  list(userBwMu=userBwMu, kernel=kern))[['mu']]
      muFun <- approxfun(Tin, xmu)
      x[['Ly']] <- lapply(1:length(x[['Ly']]), function(i)
        x[['Ly']][[i]]- muFun(x[['Lt']][[i]]))
      xmu <- muFun
    }
    
    list(x=x, mu=xmu)
  })
  
  xList <- lapply(tmp, `[[`, 'x')
  muList <- lapply(tmp, `[[`, 'mu')
  
  list(xList = xList, muList = muList)
}

## Multivariate function/scaler covariance.
# INPUTS: same as FCReg
# Output: a 4-D array containing the covariances. The first two dimensions corresponds to 
# time s and t, and the last two dimensions correspond to the variables taken covariance upon.
MvCov <- function(vars, userBwCov, outGrid, kern, measurementError=TRUE, center=TRUE, diag1D='none') {
  if (!is.list(vars) || length(vars) < 1)
    stop('`vars` needs to be a list of length >= 1')
  
  if (diag1D == 'all' && measurementError) {
    stop("Cannot assume measurement error when diag1D == 'all'")
  }
  isFuncVars <- sapply(vars, is.list)
  p <- length(isFuncVars)
  pFunc <- sum(isFuncVars)
  pScaler <- sum(!isFuncVars)
  
  if (any(isFuncVars)) {
    tAll <- do.call(c, lapply(vars[isFuncVars], function(x) unlist(x[['Lt']])))
    Tin <- sort(unique(tAll))
    
    if (missing(outGrid))
      outGrid <- Tin
    lenoutGrid <- length(outGrid)
    
  } else {
    stop('No functional observation found')
  }
  
  # First two dimensions are for s, t, and the last two dimensions are for matrix of # random variables.
  res <- array(NA, c(lenoutGrid, lenoutGrid, p, p))
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      if (j <= i) {
        use1D <- diag1D == 'all' || ( diag1D == 'cross' && j != i )
        covRes <- uniCov(vars[[i]], vars[[j]], userBwCov, outGrid, kern, 
                         rmDiag = (i == j) && measurementError, 
                         center, use1D)
        if (attr(covRes, 'covType') %in% c('FF', 'SS'))
          res[, , i, j] <- covRes
        else {
          if (nrow(covRes) == 1)   # cov(scaler, function)
            res[, , i, j] <- matrix(covRes, lenoutGrid, lenoutGrid, byrow=TRUE)
          else                     # cov(function, scaler)
            res[, , i, j] <- matrix(covRes, lenoutGrid, lenoutGrid, byrow=FALSE)
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
# center: whether to center the covariates before calculate covariance.
# use1D: whether to use 1D smoothing for estimating the diagonal covariance.
uniCov <- function(X, Y, userBwCov, outGrid, kern='gauss', rmDiag=FALSE, center=TRUE, use1D=FALSE) {
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
    Tin <- sort(unique(unlist(X[['Lt']])))
    if (center) {
      Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], Tin, Tin[1], list(userBwMu=userBwCov, kernel=kern))[['mu']]
      Ymu <- mean(Y)
    } else {
      Xmu <- rep(0, length(Tin))
      Ymu <- 0
    }
    res <- GetCrCovYZ(userBwCov, Y, Ymu, X[['Ly']], X[['Lt']], Xmu, Tin, kern)[['smoothedCC']]
    res <- as.matrix(ConvertSupport(Tin, outGrid, mu=res))
    if (flagScalerFunc) 
      res <- t(res)
    
    attr(res, 'covType') <- 'FS'
    
    # function-function  
  } else {
    TinX <- sort(unique(unlist(X[['Lt']])))
    TinY <- sort(unique(unlist(Y[['Lt']])))
    noutGrid <- length(outGrid)
    if (center) {
      if (min(TinX) > min(outGrid) || min(TinY) > min(outGrid) || 
          max(TinY) < max(outGrid) || max(TinX) < max(outGrid))
        stop('Observation time points coverage too low')
      
      Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], TinX, TinX[1],
                                  list(userBwMu=userBwCov, kernel=kern))[['mu']]
      Ymu <- GetSmoothedMeanCurve(Y[['Ly']], Y[['Lt']], TinY, TinY[1],
                                  list(userBwMu=userBwCov, kernel=kern))[['mu']]
    } else {
      Xmu <- rep(0, length(TinX))
      Ymu <- rep(0, length(TinY))
    }
    names(Xmu) <- TinX
    names(Ymu) <- TinY
    
    if (use1D) {
      Xvec <- unlist(X[['Ly']])
      Yvec <- unlist(Y[['Ly']])
      tvecX <- unlist(X[['Lt']])
      tvecY <- unlist(Y[['Lt']])
      if (!identical(tvecX, tvecY)){
        stop('Cannot use 1D covariance smoothing if the observation time points for X and Y are different')
      }
      
      ord <- order(tvecX)
      tvecX <- tvecX[ord]
      Xvec <- Xvec[ord]
      Yvec <- Yvec[ord]
      Xcent <- Xvec - Xmu[as.character(tvecX)]
      Ycent <- Yvec - Ymu[as.character(tvecX)]
      covXY <- Lwls1D(userBwCov, kern, npoly=1L, nder=0L, 
                      xin=tvecX, yin=Xcent * Ycent, 
                      win=rep(1, length(tvecX)), xout=outGrid)
      res <- matrix(NA, noutGrid, noutGrid)
      diag(res) <- covXY
    } else { # use 2D smoothing
      tmp <- GetCrCovYX(userBwCov, userBwCov, X[['Ly']], X[['Lt']], Xmu, 
                        Y[['Ly']], Y[['Lt']], Ymu, rmDiag=rmDiag, kern=kern)
      gd <- tmp[['smoothGrid']]
      res <- matrix(
              interp2lin(as.numeric(gd[, 1]), 
                         as.numeric(gd[, 2]), 
                         matrix(as.numeric(tmp[['smoothedCC']]),
                                nrow(tmp[['smoothedCC']]),
                                ncol(tmp[['smoothedCC']])), 
                         rep(as.numeric(outGrid), times=noutGrid), 
                         rep(as.numeric(outGrid), each=noutGrid)), 
              noutGrid, noutGrid)
    }
    attr(res, 'covType') <- 'FF'
  }
  
  return(res)
}

## Concurrent functional regression by imputation. This does not provide consistent estimates.
## FPCAlist: a list of functional covariates and response. Each field corresponds to a covariate. 
#            The last entry is assumed to be the response if no entry is names 'Y'.
imputeConReg <- function(FPCAlist, Z, outGrid) {
  
  if (is.null(names(FPCAlist)))
    names(FPCAlist) <- c(paste0('X', seq_len(length(FPCAlist) - 1)), 'Y')
  
  if ('Y' %in% names(FPCAlist)) {
    Yname <- 'Y'
    FPCAlist <- c(FPCAlist[names(FPCAlist) != 'Y'], FPCAlist['Y'])
  } else 
    Yname <- names(FPCAlist)[length(FPCAlist)]
  
  imputeCurves <- sapply(FPCAlist, function(x) 
    apply(fitted(x), 1, function(fit) 
      approx(x[['workGrid']], fit, outGrid)[['y']]),
    simplify='array')
  alphaBeta <- apply(imputeCurves, 1, function(XYt) {
    Yt <- XYt[, ncol(XYt)]
    designMat <- cbind(1, XYt[, -ncol(XYt), drop=FALSE], Z)
    beta_t <- qr.solve(designMat, Yt)
    return(beta_t)
  })
  beta0 <- alphaBeta[1, ]
  beta <- alphaBeta[-1, , drop=FALSE]
  
  return(list(beta0 = beta0, beta = beta, outGrid = outGrid))
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
