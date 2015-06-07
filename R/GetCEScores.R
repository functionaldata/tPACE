# obsGrid: supports fittedCov and phi. May be a truncated version.
# Assumes each tVec in t are all supported on obsGrid.
# return: ret is a 2 by n array, with the first row containing the xiEst and second row containing the xiVar.
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

  obsGrid <- round(obsGrid, 14)
  ret <- lapply(t, function(tvec) {
    ind <- match(round(tvec, 14), obsGrid)
    if (sum(is.na(ind)) != 0)
      stop('Time point not found in obsGrid.')
    
    return(list(muVec=mu[ind], phiMat=phi[ind, , drop=FALSE], Sigma_Yi=Sigma_Y[ind, ind, drop=FALSE]))
  })

  return(ret)
}


GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi,
                           verbose=FALSE) {
  if (length(yVec) == 0) {
    if (verbose)
      warnings('Empty observation found, possibly due to truncation')

    return(list(xiEst=NULL, xiVar=NULL))
  }

  Lam <- diag(lamVec)
  LamPhi <- Lam %*% t(phiMat)
  LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
  xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
  xiVar <- Lam - LamPhi %*% t(LamPhiSig)

  return(list(xiEst=xiEst, xiVar=xiVar))
}
