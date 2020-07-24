GetRho <- function(y, t, optns, mu,muwork, obsGrid, fittedCov, lambda, phi, phiwork, workgrid, sigma2) {
  
  optnsTmp <- optns
  optnsTmp$verbose <- FALSE 
  for (j in 1:2) {
    yhat <- GetCEScores(y, t, optnsTmp, mu, obsGrid, fittedCov, lambda, phi, sigma2)[3, ] 
    sigma2 <- mean(mapply(function(a, b) mean((a - b)^2, na.rm=TRUE), yhat, y), na.rm=TRUE)
  }
  ls = sapply(t,length)
  y = y[ls >0]
  t= t[ls > 0]
  
  if(optns$methodRho == 'trunc'){
    R <- sqrt((trapzRcpp(obsGrid, mu ^ 2) + sum(lambda)) / diff(range(obsGrid)))
    a1 <- -13; a2 <- -1.5
    etaCand <- seq(a1, a2, length.out=50)
    rhoCand <- exp(etaCand) * R
  }else{
    rhoCand <- seq(1,10,length.out = 50) * sigma2
  }
  
  #  cvScores <- sapply(rhoCand, cvRho, leaveOutInd=leaveOutInd, y=y, t=t, optns=optns, mu=mu, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi)
  rhoScores <- sapply(rhoCand, errorRho, y=y, t=t, optns=optns, mu=mu,muwork = muwork, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi,phiwork = phiwork,workgrid = workgrid)
  
#  plot(log(rhoCand),rhoScores,type = 'l',xlab = 'log(rho)')
  
  return(rhoCand[which.min(rhoScores)])
}


# sigma2* = max(rho, sigma2) = rho
# Get the error in total variance with 
errorRho = function(y, t, optns, mu, muwork, obsGrid, fittedCov, lambda, phi, phiwork,workgrid,rho){
  scoresObj = GetCEScores(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2=rho)
  xiEst <- t(do.call(cbind, scoresObj[1, ])) 
  K = length(lambda)
  ZMFV = xiEst[, seq_len(K), drop = FALSE] %*% t(phiwork[, seq_len(K), drop = FALSE]);   
  IM = muwork
  fittedY = t(apply( ZMFV, 1, function(x) x + IM))
  #  matplot(workgrid,t(fittedY),type="l",xlab = 'age',ylab="Fitted Trajectories")
  varY = apply(fittedY,2,var)
  intevarY = trapzRcpp(X=workgrid,Y=varY)
  return((sum(lambda)-intevarY)^2)
}

# Method of choosing rho in Bitao's paper. Get the CV score for a given rho. Correspond to getScores2.m
#cvRho <- function(rho, leaveOutInd, y, t, optns, mu, obsGrid, fittedCov, lambda, phi) {

#  Sigma_Y <- fittedCov + diag(rho, nrow(phi))

#  MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)

#  yhat <- mapply(function(yVec, muphisig, ind) 
#         GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat,
#                        muphisig$Sigma_Yi, newyInd=ind)$fittedY, 
#                y, MuPhiSig, leaveOutInd) 

#  yobs <- mapply(`[`, y, leaveOutInd)

#  return(sum((na.omit(unlist(yobs)) - unlist(yhat))^2, na.rm=TRUE))
#}


# sample one observation from each tVec in t. The 'non-random' sampling is for testing against Matlab.
RandTime <- function(t, isRandom=TRUE) {
  ni <- sapply(t, length)
  if (all(ni <= 1)) 
    stop('None of the individuals have >= 2 observations. Cannot use rho')
  
  if (isRandom) {
    ind <- sapply(ni, sample, size=1)
  } else {
    ind <- ((1000 + 1:length(t)) %% ni) + 1
  }
  
  return(ind)
}
