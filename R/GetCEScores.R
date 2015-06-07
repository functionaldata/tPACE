GetCEScores <- function(y, t, optns, mu, obsGrid, lambda, phi, sigma2) {

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

  ret <- lapply(t, function(tvec) {
    ind <- match(tvec, obsGrid)
    if (sum(is.na(ind)) != 0)
      stop('Time point not found in obsGrid.')
    
    return(list(muVec=mu[ind], phiMat=phi[ind, ], Sigma_Yi=Sigma_Y[ind, ind]))
  })

  return(ret)
}


GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi) {

  Lam <- diag(lamVec)
  LamPhi <- Lam %*% t(phiMat)
  LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
  xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
  xiVar <- Lam - LamPhi %*% t(LamPhiSig)

  return(xiEst=xiEst, xiVar=xiVar)
}
