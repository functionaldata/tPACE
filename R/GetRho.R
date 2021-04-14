GetRho <- function(y, t, optns, mu, muwork, obsGrid, fittedCov, lambda, phi, phiwork, workgrid, sigma2,idxRho) {
  
  optnsTmp <- optns
  optnsTmp$verbose <- FALSE
  
  if((optns$methodRho != 'ridge')&& (optns$methodRho != 'RidgeLeaveOneOutPerSubjectCV')){ #If methodRho is ridge then this is not used so we can speed up the code
    for (j in 1:2) {
      yhat <- GetCEScores(y, t, optnsTmp, mu, obsGrid, fittedCov, lambda, phi, sigma2)[3, ] 
      sigma2 <- mean(mapply(function(a, b) mean((a - b)^2, na.rm=TRUE), yhat, y), na.rm=TRUE)
    }
  }
  
  ls = sapply(t,length) 
  y = y[ls >0]
  t= t[ls > 0]
  
  if(optns$methodRho == 'ridge'){
    n=length(y)
    rhoOut=optimize(errorRhoCV,lower=10^(-9),upper=10*sigma2,y=y[idxRho], t=t[idxRho], optns=optns, mu=mu,muwork = muwork, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi,phiwork = phiwork,workgrid = workgrid)$minimum
  }else if(optns$methodRho == 'RidgeLeaveOneOutPerSubjectCV'){
    rhoOut <- optimize(errorRhoLeaveOneOutPerSubjectCV,lower=10^(-9),upper=(10*sigma2),y=y[idxRho], t=t[idxRho], optns=optns, mu=mu,muwork = muwork, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi,phiwork = phiwork,workgrid = workgrid)$minimum
  }else{
    rhoCand <- seq(1,10,length.out = 50) * sigma2
    rhoScores <- sapply(rhoCand, errorRho, y=y[idxRho], t=t[idxRho], optns=optns, mu=mu,muwork = muwork, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi,phiwork = phiwork,workgrid = workgrid)
    rhoOut=rhoCand[which.min(rhoScores)]
  }
  return(rhoOut)
}

errorRhoLeaveOneOutPerSubjectCV = function(y, t, optns, mu, muwork, obsGrid, fittedCov, lambda, phi, phiwork,workgrid,rho){
  rep = 1 # rep-fold cross validation.
  err <- rep(0, rep)
  for(i in c(1:rep)){
    set.seed(i+2020)
    TestIndex <- RandTime(t, isRandom=TRUE)
    yTest <- lapply(c(1:length(t)), function(j){return(y[[j]][TestIndex[j]])})
    tTest <- lapply(c(1:length(t)), function(j){return(t[[j]][TestIndex[j]])})
    yTrain <- lapply(c(1:length(t)), function(j){return(y[[j]][-TestIndex[j]])})
    tTrain <- lapply(c(1:length(t)), function(j){return(t[[j]][-TestIndex[j]])})
    
    scoresObj = GetCEScores(yTrain, tTrain, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2=rho)
    xiEst <- t(do.call(cbind, scoresObj[1, ]))
    K = length(lambda)
    xiEst=xiEst[, seq_len(K), drop = FALSE]
    
    for(j in 1:length(t)){
      MuObs = ConvertSupport(fromGrid = workgrid, toGrid = tTest[[j]], mu = muwork)
      PhiObs = ConvertSupport(fromGrid = workgrid, toGrid = tTest[[j]], phi = phiwork)
      if(length(tTest[[j]])==1){
        PhiObs = PhiObs[seq_len(K), drop = FALSE]
        ZMFV =xiEst[j,] %*% PhiObs
      }else{
        PhiObs = PhiObs[,seq_len(K), drop = FALSE]
        ZMFV =xiEst[j,] %*% t(PhiObs)
      }
      IM = MuObs
      fittedY = t(apply( ZMFV, 1, function(x) x + IM))
      err[i]=err[i]+sum(fittedY-yTest[[j]])^2/length(y)
    }
  }
  return(sum(err)/rep)
}



errorRhoCV = function(y, t, optns, mu, muwork, obsGrid, fittedCov, lambda, phi, phiwork,workgrid,rho){
  predictionError=0
  scoresObj = GetCEScores(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2=rho)
  xiEst <- t(do.call(cbind, scoresObj[1, ]))
  K = length(lambda)
  xiEst=xiEst[, seq_len(K), drop = FALSE]

  for(j in 1:length(t)){
    MuObs = ConvertSupport(fromGrid = workgrid, toGrid = t[[j]], mu = muwork)
    PhiObs = ConvertSupport(fromGrid = workgrid, toGrid = t[[j]], phi = phiwork)
    if(length(t[[j]])==1){
      PhiObs = PhiObs[seq_len(K), drop = FALSE]
      ZMFV =xiEst[j,] %*% PhiObs
    }else{
      PhiObs = PhiObs[,seq_len(K), drop = FALSE]
      ZMFV =xiEst[j,] %*% t(PhiObs)
    }
    IM = MuObs
    fittedY = t(apply( ZMFV, 1, function(x) x + IM))
    predictionError=predictionError+sum(fittedY-y[[j]])^2/length(y)
  }
  return(predictionError)
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
