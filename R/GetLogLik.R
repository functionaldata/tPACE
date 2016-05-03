# k: input denoting number of components used
# returns -2 times log-likelihood
GetLogLik = function(fpcaObj, k, Ly = NULL, Lt = NULL){
  if(fpcaObj$optns$lean == TRUE && (is.null(Ly) || is.null(Lt))){
    stop("Option lean is TRUE, need input data Ly and measurement time list Lt to calculate log-likelihood.")
  }
  if(fpcaObj$optns$lean == FALSE){ # when input data is in fpcaObj
    Ly <- fpcaObj$inputData$Ly
    Lt <- fpcaObj$inputData$Lt
  }
  lambda = fpcaObj$lambda[1:k]
  sigma2 = fpcaObj$sigma2
  if(is.null(sigma2) && fpcaObj$optns$dataType == "Dense"){
    ymat = matrix(unlist(Ly),nrow=length(Ly), byrow=TRUE)
    sddiag = sqrt(diag(var(ymat)))
    sigma2 = sddiag*1e-4
    sigma2 = ConvertSupport(fromGrid = fpcaObj$obsGrid, toGrid = fpcaObj$workGrid, mu = sigma2)
  }
  logLik = 0
  phi = fpcaObj$phi[,1:k, drop=FALSE]

  if(fpcaObj$optns$dataType %in% c('Dense'
    #, 'DenseWithMV' # need extra imputation step
    )){
  	if(k == 1){
  	  Sigma_y = phi %*% (lambda*diag(k)) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
  	} else {
      Sigma_y = phi %*% diag(lambda) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
    }
    detSigma_y = prod(c(lambda,rep(0,nrow(phi)-k))[1:length(lambda)]+sigma2)
    #detSigma_y = det(Sigma_y)
    if(detSigma_y == 0){
      logLik = NULL
      return(logLik)
    }
    # calculate loglikelihood via matrix multiplication
    ymatcenter = matrix(unlist(Ly)-fpcaObj$mu, nrow = length(Ly), byrow = TRUE)
    svd_Sigma_y = svd(Sigma_y)
    Sigma_y_inv = svd_Sigma_y$v %*% diag(1/svd_Sigma_y$d) %*% t(svd_Sigma_y$u)
  	logLik = sum(diag(t(Sigma_y_inv %*% t(ymatcenter)) %*% t(ymatcenter))) + length(Ly)*log(detSigma_y)
    return(logLik)
  } else { # Sparse case
    if(is.null(sigma2)){ sigma2 <- fpcaObj$rho }
    if(fpcaObj$optns$error == TRUE && sigma2 <= fpcaObj$rho){
      # especially for the case when sigma2 is estimated to be <=0 and set to 1e-6
      sigma2 <- fpcaObj$rho
    }
    for(i in 1:length(Ly)){
      if(length(Lt[[i]]) == 1){
        phi_i = t(as.matrix(ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
                               phi = phi)))        
      } else {
        phi_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
                               phi = phi)
      }
      mu_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
        mu = fpcaObj$mu)
      if(k == 1){
        Sigma_yi = phi_i %*% (lambda*diag(k)) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      } else{
        Sigma_yi = phi_i %*% diag(lambda) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      }
      detSigma_yi = det(Sigma_yi)
      if(detSigma_yi == 0){
        logLik = NULL
        return(logLik)
      }
      invtempi = solve(Sigma_yi, Ly[[i]] - mu_i)
      logLik = logLik + log(detSigma_yi) + invtempi %*% (Ly[[i]] - mu_i)
    }
    return(logLik)
  }
}
