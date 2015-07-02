# obsGrid: supports fittedCov and phi. May be a truncated version.
# Assumes each tVec in t are all supported on obsGrid.
# return: ret is a 3 by n array, with the first row containing the xiEst, second row containing the xiVar, and the third containing the fitted values. 
GetCEScores <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {

  if (length(lambda) != ncol(phi))
    stop('No of eigenvalues is not the same as the no of eigenfunctions.')

  if (is.null(sigma2))
    sigma2 <- 0

  Sigma_Y <- fittedCov + diag(sigma2, nrow(phi))

  MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)
  ret <- mapply(function(yVec, muphisig) 
         GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat, muphisig$Sigma_Yi), 
         y, MuPhiSig) 

  return(ret)
}


GetMuPhiSig <- function(t, obsGrid, mu, phi, Sigma_Y) {

  obsGrid <- signif(obsGrid, 14)
  ret <- lapply(t, function(tvec) {
    ind <- match(signif(tvec, 14), obsGrid)
    if (sum(is.na(ind)) != 0)
      stop('Time point not found in obsGrid.')
    
    return(list(muVec=mu[ind], phiMat=phi[ind, , drop=FALSE], Sigma_Yi=Sigma_Y[ind, ind, drop=FALSE]))
  })

  return(ret)
}


GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd=NULL,
                            verbose=FALSE) {

  if (length(yVec) == 0) {
    if (verbose)
      warnings('Empty observation found, possibly due to truncation')

      return(list(xiEst=NULL, xiVar=NULL, fittedY=NULL))
  }

# When an individual has only one observation, the leave-one-out predicted Y is NA.
  if (length(yVec) == 1 && !is.null(newyInd)) {
    newPhi <- matrix(NA, ncol=length(lamVec))
    newMu <- NA
  }

  if (!is.null(newyInd) && length(yVec) != 1) {
    newy <- yVec[newyInd]
    newPhi <- phiMat[newyInd, , drop=FALSE]
    newMu <- muVec[newyInd]

    yVec <- yVec[-newyInd]
    muVec <- muVec[-newyInd]
    phiMat <- phiMat[-newyInd, , drop=FALSE]
    Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop=FALSE]
  }

  Lam <- diag(lamVec)
  LamPhi <- Lam %*% t(phiMat)
  LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
  xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
  xiVar <- Lam - LamPhi %*% t(LamPhiSig)


  fittedY <- if(is.null(newyInd)) 
    muVec + phiMat %*% xiEst else 
      newMu + newPhi %*% xiEst

  ret <- list(xiEst=xiEst, xiVar=xiVar, fittedY=fittedY)

  return(ret)
}
