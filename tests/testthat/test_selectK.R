library(testthat)
library(mvtnorm)
devtools::load_all()

test_that("SelectK works for dense FPCA Wiener process with measurement error",{
  n = 100; K = 10; pts = seq(0, 1, length=50)
  lambda <- (1 / (1:K - 1/2) / pi)^2
  TrueFVE <- cumsum(lambda)/sum(lambda) * 100
  phi <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
  set.seed(1)
  scores <- t(diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
  sd <- 0.01
  # generate Wiener process
  test1 <- t(phi %*% t(scores)) # Dense data
  test2 <- test1 + matrix(rnorm(nrow(test1)*ncol(test1),mean=0,sd=sd),nrow=nrow(test1))
  test4 <- Sparsify(test2, pts, rep(length(pts),n)) # Dense Data with measurement errors
  Trueloglik <- rep(0, K)
  for(i in 1:K){
    # marginal likelihood
    Sigma_y <- phi[,1:i,drop=F] %*% (lambda[1:i]*diag(i)) %*% t(phi[,1:i,drop=F]) + sd^2 * diag(nrow(phi[,1:i,drop=F]))
    detSigma_y <- det(Sigma_y)
    for(j in 1:n){
      invtempsub = solve(Sigma_y, test4$yList[[j]] - 0)
      Trueloglik[i] <- Trueloglik[i] + dmvnorm(x=test4$yList[[j]],mean=rep(0,nrow(phi)), 
                                               sigma=Sigma_y, log=TRUE)
    }
    # conditional likelihood
    #Trueloglik[i] <- -(sum(diag((test2 - t(phi[,1:i] %*% t(scores[,1:i]))) %*% t(test2 - t(phi[,1:i] %*% t(scores[,1:i])))))/(2*sd^2) + 
    #                     length(pts)/2*log(2*pi) + length(pts)/2*log(sd^2))
  }
  TrueAIC <- 2*(1:K) - 2*Trueloglik 
  TrueBIC <- log(n)*(1:K) - 2*Trueloglik
  # FPCA
  #optnDenseError <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
  optnDenseErrorAIC <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'AIC', outPercent=c(0, 1), verbose=TRUE))
  optnDenseErrorBIC <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'BIC', outPercent=c(0, 1), verbose=TRUE))  
  optnDenseErrorFVE <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'FVE', FVEthreshold = 0.85, outPercent=c(0, 1), verbose=TRUE, lean = FALSE))  
  test4FPCA_AIC <- FPCA(test4$yList, test4$tList, optnDenseErrorAIC)
  test4FPCA_BIC <- FPCA(test4$yList, test4$tList, optnDenseErrorBIC)
  test4FPCA_FVE <- FPCA(test4$yList, test4$tList, optnDenseErrorFVE)
  expect_equal(test4FPCA_FVE$selectK, min(which(TrueFVE >= 85)))
  expect_equal(test4FPCA_AIC$selectK, which(TrueAIC == min(TrueAIC)))
  expect_equal(test4FPCA_BIC$selectK, which(TrueBIC == min(TrueBIC)))
})

#test_that("SelectK works for sparse FPCA Wiener process with measurement error",{
#  n = 1000; K = 5; pts = seq(0, 1, length=50)
#  lambda <- (1 / (1:K - 1/2) / pi)^2
#  TrueFVE <- cumsum(lambda)/sum(lambda) * 100
#  phi <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
#  set.seed(1)
#  scores <- t(diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
#  sd <- 0.01
#  # generate Wiener process
#  test1 <- t(phi %*% t(scores)) # Dense data
#  test2 <- test1 + matrix(rnorm(nrow(test1)*ncol(test1),mean=0,sd=sd),nrow=nrow(test1))
#  test3 <- Sparsify(test2, pts, 4:10) # Sparse Data with measurement errors
#  Trueloglik <- rep(0, K)
#  for(i in 1:K){
#    # marginal likelihood
#    for(j in 1:n){
#      phi_i = ConvertSupport(fromGrid = pts, toGrid = test3$tList[[j]], phi = phi)
#      Sigma_yi <- phi_i[,1:i,drop=F] %*% (lambda[1:i]*diag(i)) %*% t(phi_i[,1:i,drop=F]) + sd^2 * diag(nrow(phi_i[,1:i,drop=F]))
#      detSigma_yi <- det(Sigma_yi)
#      invtempsub = solve(Sigma_yi, test3$yList[[j]] - 0)
#      Trueloglik[i] <- Trueloglik[i] + dmvnorm(x=test3$yList[[j]],mean=rep(0,nrow(phi_i)), 
#                                               sigma=Sigma_yi, log=TRUE)
#    }
#    # conditional likelihood
#    #Trueloglik[i] <- -(sum(diag((test2 - t(phi[,1:i] %*% t(scores[,1:i]))) %*% 
#    #                     t(test2 - t(phi[,1:i] %*% t(scores[,1:i])))))/(2*sd^2) - 
#    #                     n*length(pts)/2*log(2*pi) - n*length(pts)/2*log(sd^2))
#  }
#  TrueAIC <- 2*(1:K) - 2*Trueloglik 
#  TrueBIC <- log(n)*(1:K) - 2*Trueloglik
#  # FPCA
#  #optnSparseError <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=FALSE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
#  optnSparseError_AIC <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='AIC', outPercent=c(0, 1), verbose=TRUE))
#  optnSparseError_BIC <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='BIC', outPercent=c(0, 1), verbose=TRUE))
#  optnSparseError_FVE <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='FVE', FVEthreshold = 0.85, outPercent=c(0, 1), verbose=TRUE))
#  test3FPCA_AIC <- FPCA(test3$yList, test3$tList, optnSparseError_AIC)
#  test3FPCA_BIC <- FPCA(test3$yList, test3$tList, optnSparseError_BIC)
#  test3FPCA_FVE <- FPCA(test3$yList, test3$tList, optnSparseError_FVE)
#  expect_equal(test3FPCA_AIC$selectK, which(TrueAIC == min(TrueAIC)))
#  expect_equal(test3FPCA_BIC$selectK, which(TrueBIC == min(TrueBIC)))
#  expect_equal(test3FPCA_FVE$selectK, min(which(TrueFVE >= 85)))
#})

test_that("SelectK works outside FPCA function for Dense data and gives the same output with equivalent options",{
  n = 100; K = 10; pts = seq(0, 1, length=50)
  lambda <- (1 / (1:K - 1/2) / pi)^2
  TrueFVE <- cumsum(lambda)/sum(lambda) * 100
  phi <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
  set.seed(1)
  scores <- t(diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
  sd <- 0.01
  # generate Wiener process
  test1 <- t(phi %*% t(scores)) # Dense data
  test2 <- test1 + matrix(rnorm(nrow(test1)*ncol(test1),mean=0,sd=sd),nrow=nrow(test1))
  test4 <- Sparsify(test2, pts, rep(length(pts),n)) # Dense Data with measurement errors
  Trueloglik <- rep(0, K)
  for(i in 1:K){
    # marginal likelihood
    Sigma_y <- phi[,1:i,drop=F] %*% (lambda[1:i]*diag(i)) %*% t(phi[,1:i,drop=F]) + sd^2 * diag(nrow(phi[,1:i,drop=F]))
    detSigma_y <- det(Sigma_y)
    for(j in 1:n){
      invtempsub = solve(Sigma_y, test4$yList[[j]] - 0)
      Trueloglik[i] <- Trueloglik[i] + dmvnorm(x=test4$yList[[j]],mean=rep(0,nrow(phi)), 
                                               sigma=Sigma_y, log=TRUE)
    }
    # conditional likelihood
    #Trueloglik[i] <- -(sum(diag((test2 - t(phi[,1:i] %*% t(scores[,1:i]))) %*% t(test2 - t(phi[,1:i] %*% t(scores[,1:i])))))/(2*sd^2) + 
    #                     length(pts)/2*log(2*pi) + length(pts)/2*log(sd^2))
  }
  TrueAIC <- 2*(1:K) - 2*Trueloglik 
  TrueBIC <- log(n)*(1:K) - 2*Trueloglik
  # FPCA
  #optnDenseError <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
  optnDenseErrorAIC <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'AIC', outPercent=c(0, 1), verbose=TRUE))
  optnDenseErrorBIC <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'BIC', outPercent=c(0, 1), verbose=TRUE))  
  optnDenseErrorFVE <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = 'FVE', FVEthreshold = 0.85, outPercent=c(0, 1), verbose=TRUE, lean = FALSE))  
  test4FPCA_AIC <- FPCA(test4$yList, test4$tList, optnDenseErrorAIC)
  test4FPCA_BIC <- FPCA(test4$yList, test4$tList, optnDenseErrorBIC)
  test4FPCA_FVE <- FPCA(test4$yList, test4$tList, optnDenseErrorFVE)
  # the following default FPCA run does not select number of components
  optnDenseError <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE, lean = FALSE))
  test4FPCA <- FPCA(test4$yList, test4$tList, optnDenseError)
  outAIC <- SelectK(fpcaObj = test4FPCA, criterion = 'AIC')$k
  outBIC <- SelectK(fpcaObj = test4FPCA, criterion = 'BIC')$k
  outFVE <- SelectK(fpcaObj = test4FPCA, criterion = 'FVE', FVEthreshold = 0.85)$k
  expect_equal(test4FPCA_AIC$selectK, outAIC)
  expect_equal(test4FPCA_BIC$selectK, outBIC)
  expect_equal(test4FPCA_FVE$selectK, outFVE)
  # the option fixedK is tested below
  FixK <- 3
  optnDenseErrorfixedK <- SetOptions(test4$yList, test4$tList, list(dataType='Dense', error=TRUE, kernel='epan', selectionMethod = FixK, outPercent=c(0, 1), verbose=TRUE, lean = FALSE))  
  test4FPCA_fixedK <- FPCA(test4$yList, test4$tList, optnDenseErrorfixedK)
  expect_equal(test4FPCA_fixedK$selectK, FixK)
  # we can get the same results if using SelectK outside FPCA
  outfixedK <- SelectK(fpcaObj = test4FPCA, criterion = FixK)$k
  expect_equal(outfixedK, FixK)
})

test_that("SelectK works outside FPCA function for Sparse data and gives the same output with equivalent options",{
  n = 1000; K = 5; pts = seq(0, 1, length=50)
  lambda <- (1 / (1:K - 1/2) / pi)^2
  TrueFVE <- cumsum(lambda)/sum(lambda) * 100
  phi <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
  set.seed(1)
  scores <- t(diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
  sd <- 0.01
  # generate Wiener process
  test1 <- t(phi %*% t(scores)) # Dense data
  test2 <- test1 + matrix(rnorm(nrow(test1)*ncol(test1),mean=0,sd=sd),nrow=nrow(test1))
  test3 <- Sparsify(test2, pts, 4:10) # Sparse Data with measurement errors
  Trueloglik <- rep(0, K)
  for(i in 1:K){
    # marginal likelihood
    for(j in 1:n){
      phi_i = ConvertSupport(fromGrid = pts, toGrid = test3$tList[[j]], phi = phi)
      Sigma_yi <- phi_i[,1:i,drop=F] %*% (lambda[1:i]*diag(i)) %*% t(phi_i[,1:i,drop=F]) + sd^2 * diag(nrow(phi_i[,1:i,drop=F]))
      detSigma_yi <- det(Sigma_yi)
      invtempsub = solve(Sigma_yi, test3$yList[[j]] - 0)
      Trueloglik[i] <- Trueloglik[i] + dmvnorm(x=test3$yList[[j]],mean=rep(0,nrow(phi_i)), 
                                               sigma=Sigma_yi, log=TRUE)
    }
    # conditional likelihood
    #Trueloglik[i] <- -(sum(diag((test2 - t(phi[,1:i] %*% t(scores[,1:i]))) %*% 
    #                     t(test2 - t(phi[,1:i] %*% t(scores[,1:i])))))/(2*sd^2) - 
    #                     n*length(pts)/2*log(2*pi) - n*length(pts)/2*log(sd^2))
  }
  TrueAIC <- 2*(1:K) - 2*Trueloglik 
  TrueBIC <- log(n)*(1:K) - 2*Trueloglik
  # FPCA
  #optnSparseError <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=FALSE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
  optnSparseError_AIC <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='AIC', outPercent=c(0, 1), verbose=TRUE))
  optnSparseError_BIC <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='BIC', outPercent=c(0, 1), verbose=TRUE))
  optnSparseError_FVE <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', selectionMethod='FVE', FVEthreshold = 0.85, outPercent=c(0, 1), verbose=TRUE))
  test3FPCA_AIC <- FPCA(test3$yList, test3$tList, optnSparseError_AIC)
  test3FPCA_BIC <- FPCA(test3$yList, test3$tList, optnSparseError_BIC)
  test3FPCA_FVE <- FPCA(test3$yList, test3$tList, optnSparseError_FVE)
  # the following default FPCA run does not select number of components
  optnSparseError <- SetOptions(test3$yList, test3$tList, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE, lean = FALSE))
  test3FPCA <- FPCA(test3$yList, test3$tList, optnSparseError)
  outAIC <- SelectK(fpcaObj = test3FPCA, criterion = 'AIC')$k
  outBIC <- SelectK(fpcaObj = test3FPCA, criterion = 'BIC')$k
  outFVE <- SelectK(fpcaObj = test3FPCA, criterion = 'FVE', FVEthreshold = 0.85)$k
  expect_equal(test3FPCA_AIC$selectK, outAIC)
  expect_equal(test3FPCA_BIC$selectK, outBIC)
  expect_equal(test3FPCA_FVE$selectK, outFVE)
})
