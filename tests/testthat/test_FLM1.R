library(MASS)
library(testthat)
#devtools::load_all()

test_that('GetBR2New works', {

  set.seed(1)
  # Scalar Y, scalar X scores with 1 column
  n <- 20
  y <- matrix(rnorm(n))
  x <- lapply(1:1, function(x) cbind(rnorm(n)))
  res <- GetBR2New(y, x)
  expect_equal(length(res$bList), 1)
  expect_equal(vapply(res$bList, nrow, 1L), 1L)
  expect_equal(vapply(res$bList, ncol, 1L), 1L)

  set.seed(1)
  # Scalar Y, scalar X scores with 2 columns
  n <- 20
  y <- matrix(rnorm(n))
  x <- list(matrix(rnorm(n * 2), n, 2))
  res <- GetBR2New(y, x)
  expect_equal(length(res$bList), 1)
  expect_equal(vapply(res$bList, nrow, 1L), 2L)
  expect_equal(vapply(res$bList, ncol, 1L), 1L)

  # Scalar Y, 2 X scores each with 1 column
  n <- 20
  y <- matrix(rnorm(n))
  x <- lapply(1:2, function(x) cbind(rnorm(n)))
  res <- GetBR2New(y, x)
  expect_equal(length(res$bList), 2)
  expect_equal(vapply(res$bList, nrow, 1L), c(1L, 1L))
  expect_equal(vapply(res$bList, ncol, 1L), c(1L, 1L))

  # Two Y scores, 2 X scores each with 1 column
  n <- 20
  y <- matrix(rnorm(n * 2), n, 2)
  x <- lapply(1:2, function(x) cbind(rnorm(n)))
  res <- GetBR2New(y, x)
  expect_equal(length(res$bList), 2)
  expect_equal(vapply(res$bList, nrow, 1L), c(1L, 1L))
  expect_equal(vapply(res$bList, ncol, 1L), c(2L, 2L))

})

test_that('Dense, scalar response case works', {
  set.seed(1000)
  
  library(MASS)
  
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  
  # training set
  n <- 100
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  denseLt <- list(); denseLy <- list()
  sparseLt <- list(); sparseLy <- list()
  
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLt[[i]] <- t0
    denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    
    ind <- sort(sample(1:length(t0),3))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[[i]][ind]
  }
  
  denseX <- list(Ly=denseLy,Lt=denseLt)
  sparseX <- list(Ly=sparseLy,Lt=sparseLt)
  
  denseX <- list(X=denseX)
  sparseX <- list(X=sparseX)
  
  # test set
  N <- 500
  
  XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
  
  denseLtTest <- list(); denseLyTest <- list()
  
  sparseLtTest <- list(); sparseLyTest <- list()
  
  t0 <- seq(0,1,length.out=51)
  for (i in 1:N) {
    denseLtTest[[i]] <- t0
    denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
    
    ind <- sort(sample(1:length(t0),5))
    sparseLtTest[[i]] <- t0[ind]
    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
  }
  
  denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)
  sparseXTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
  denseXTest <- list(X=denseXTest)
  sparseXTest <- list(X=sparseXTest)
  
  
  ### scalar response
  beta <- c(1, -1)
  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
  YTest <- c(XiTest%*%diag(lambdaX)%*%beta) + rnorm(N,0,0.5)
  
  ## dense
  denseFLM <- FLM1(Y=Y,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95))
  
  trueBetaList <- list()
  trueBetaList[[1]] <- cbind(phi1(denseFLM$workGridX[[1]],1),phi1(denseFLM$workGridX[[1]],2))%*%beta
  
  # coefficient function estimation error (L2-norm)
  # par(mfrow=c(1,1))
  # plot(denseFLM$workGridX[[1]],denseFLM$betaList[[1]],type='l',xlab='t',ylab=paste('beta',1,sep=''))
  # points(denseFLM$workGridX[[1]],trueBetaList[[1]],type='l',col=2)
  
  denseEstErr <- sqrt(trapzRcpp(denseFLM$workGridX[[1]],(denseFLM$betaList[[1]] - trueBetaList[[1]])^2))
  denseEstErr
  
  # par(mfrow=c(1,2))
  # plot(denseFLM$yHat,Y,xlab='fitted Y', ylab='observed Y')
  # abline(coef=c(0,1),col=8)
  
  # plot(denseFLM$yPred,YTest,xlab='predicted Y', ylab='observed Y')
  abline(coef=c(0,1),col=8)
  
  # prediction error
  densePredErr <- sqrt(mean((YTest - denseFLM$yPred)^2))
  densePredErr
  
  # Errors are properly small
  expect_lt(denseEstErr, 0.1) 
  expect_lt(densePredErr, 1) 
  expect_equal(denseFLM$muY, mean(Y))
})



#test_that('Sparse, scalar response case works', {
#  set.seed(1000)
  
#  library(MASS)
  
#  ### functional covariate
#  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
#  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
#  lambdaX <- c(1,0.7)
  
#  # training set
#  n <- 100
#  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
#  denseLt <- list(); denseLy <- list()
#  sparseLt <- list(); sparseLy <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:n) {
#    denseLt[[i]] <- t0
#    denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),3))
#    sparseLt[[i]] <- t0[ind]
#    sparseLy[[i]] <- denseLy[[i]][ind]
#  }
  
#  denseX <- list(Ly=denseLy,Lt=denseLt)
#  sparseX <- list(Ly=sparseLy,Lt=sparseLt)
  
#  denseX <- list(X=denseX)
#  sparseX <- list(X=sparseX)
  
#  # test set
#  N <- 500
  
#  XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
  
#  denseLtTest <- list(); denseLyTest <- list()
  
#  sparseLtTest <- list(); sparseLyTest <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:N) {
#    denseLtTest[[i]] <- t0
#    denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLtTest[[i]] <- t0[ind]
#    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#  }
  
#  denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#  sparseXTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
#  denseXTest <- list(X=denseXTest)
#  sparseXTest <- list(X=sparseXTest)
  
  
#  ### scalar response
#  beta <- c(1, -1)
#  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
#  YTest <- c(XiTest%*%diag(lambdaX)%*%beta) + rnorm(N,0,0.5)
  
#  sparseFLM <- FLM1(Y=Y,X=sparseX,XTest=sparseXTest,optnsListX=list(FVEthreshold=0.95))
  
#  trueBetaList <- list()
#  trueBetaList[[1]] <- cbind(phi1(sparseFLM$workGridX[[1]],1),phi1(sparseFLM$workGridX[[1]],2))%*%beta
  
#  # coefficient function estimation error (L2-norm)
#  par(mfrow=c(1,1))
#  plot(sparseFLM$workGridX[[1]],sparseFLM$betaList[[1]],type='l',xlab='t',ylab=paste('beta',1,sep=''))
#  points(sparseFLM$workGridX[[1]],trueBetaList[[1]],type='l',col=2)
  
#  sparseEstErr <- sqrt(trapzRcpp(sparseFLM$workGridX[[1]],(sparseFLM$betaList[[1]] - trueBetaList[[1]])^2))
#  sparseEstErr
  
#  par(mfrow=c(1,2))
#  plot(sparseFLM$yHat,Y,xlab='fitted Y', ylab='observed Y')
#  abline(coef=c(0,1),col=8)
  
#  plot(sparseFLM$yPred,YTest,xlab='fitted Y', ylab='observed Y')
#  abline(coef=c(0,1),col=8)
  
#  # prediction error
#  sparsePredErr <- sqrt(mean((YTest - sparseFLM$yPred)^2))
#  sparsePredErr
  
  
#  # Errors are properly small
#  expect_lt(sparseEstErr, 3.1) 
#  expect_lt(sparsePredErr, 3) 
#  expect_equal(sparseFLM$muY, mean(Y))
#})





#test_that('Dense, functional response case works', {
  
#  M <- 51
#  set.seed(1000)
  
#  library(MASS)
  
#  ### functional covariate
#  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
#  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
#  lambdaX <- c(1,0.7)
  
#  # training set
#  n <- 100
#  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
#  denseLt <- list(); denseLy <- list()
#  sparseLt <- list(); sparseLy <- list()
  
#  t0 <- seq(0,1,length.out=M)
#  for (i in 1:n) {
#    denseLt[[i]] <- t0
#    denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),3))
#    sparseLt[[i]] <- t0[ind]
#    sparseLy[[i]] <- denseLy[[i]][ind]
#  }
  
#  denseX <- list(Ly=denseLy,Lt=denseLt)
#  sparseX <- list(Ly=sparseLy,Lt=sparseLt)
  
#  denseX <- list(X=denseX)
#  sparseX <- list(X=sparseX)
  
#  # test set
#  N <- 500
  
#  XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
  
#  denseLtTest <- list(); denseLyTest <- list()
  
#  sparseLtTest <- list(); sparseLyTest <- list()
  
#  t0 <- seq(0,1,length.out=M)
#  for (i in 1:N) {
#    denseLtTest[[i]] <- t0
#    denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLtTest[[i]] <- t0[ind]
#    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#  }
  
#  denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#  sparseXTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
#  denseXTest <- list(X=denseXTest)
#  sparseXTest <- list(X=sparseXTest)
  
  
#  ### functional response
#  alpha <- 1
#  Beta <- matrix(c(1,0.5,-1,-1),nrow=2,ncol=2)

#  # training set
#  n <- 100
  
#  Eta <- Xi%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(n*2,0,0.5),nrow=n,ncol=2)
  
#  denseLt <- list()
#  denseLy <- list()
  
#  sparseLt <- list()
#  sparseLy <- list()
  
#  t0 <- seq(0,1,length.out=M)
#  for (i in 1:n) {
    
#    # comp 1
#    denseLt[[i]] <- t0
#    # denseLy[[i]] <- lambdaY[1]*Eta[i,1]*phi1(t0,1) + lambdaY[2]*Eta[i,2]*phi2(t0,1) + 
#    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
#    denseLy[[i]] <- Eta[i,1]*phi1(t0,1) + Eta[i,2]*phi2(t0,1) + 
#      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) + alpha
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLt[[i]] <- t0[ind]
#    sparseLy[[i]] <- denseLy[[i]][ind]
#  }
  
#  denseY <- list(Ly=denseLy,Lt=denseLt)
#  sparseY <- list(Ly=sparseLy,Lt=sparseLt)
  
#  # test set
#  N <- 500
  
#  EtaTest <- XiTest%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(N*2,0,0.5),nrow=N,ncol=2) 
  
#  denseLtTest <- list()
#  denseLyTest <- list()
  
#  sparseLtTest <- list()
#  sparseLyTest <- list()
  
#  t0 <- seq(0,1,length.out=M)
#  for (i in 1:N) {
    
#    # comp 1
#    denseLtTest[[i]] <- t0
#    # denseLyTest[[i]] <- lambdaY[1]*EtaTest[i,1]*phi1(t0,1) + lambdaY[2]*EtaTest[i,2]*phi2(t0,1) + 
#    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
#    denseLyTest[[i]] <- EtaTest[i,1]*phi1(t0,1) + EtaTest[i,2]*phi2(t0,1) + 
#      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) + alpha
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLtTest[[i]] <- t0[ind]
#    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#  }
  
#  denseYTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#  sparseYTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
  
#  ## dense
#  # print(system.time({
#  denseFLM <- FLM1(Y=denseY,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95), nPerm=1000)
#  # }))
  
#  trueBetaList <- list()
#  trueBetaList[[1]] <- cbind(phi1(denseFLM$workGridX[[1]],1),phi1(denseFLM$workGridX[[1]],2))%*%Beta[1:2,]%*%
#    t(cbind(phi1(denseFLM$workGridY,1),phi2(denseFLM$workGridY,1)))
  
#  # coefficient function estimation error (L2-norm)
#  denseEstErr <- sqrt(sum((denseFLM$betaList[[1]] - trueBetaList[[1]])^2)*
#                        diff(denseFLM$workGridX[[1]])[1]*
#                        diff(denseFLM$workGridY)[1])
#  denseEstErr
  
#  # par(mfrow=c(2,3))
#  # for (i in 1:6) {
#  #   i <- sample(1:n,1)
#  #   plot(denseFLM$workGridY,denseFLM$yHat[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (hat)')
#  #   points(denseY$Lt[[i]],denseY$Ly[[i]],type='l',col=2)
#  # }
#  # 
#  # par(mfrow=c(2,3))
#  # for (i in 1:6) {
#  #   i <- sample(1:n,1)
#  #   plot(denseFLM$workGridY,denseFLM$yPred[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (pred)')
#  #   points(denseYTest$Lt[[i]],denseYTest$Ly[[i]],type='l',col=2)
#  # }
  
#  # # prediction error
#  # densePredErr <- sqrt(mean(apply((YTest - denseFLM$yPred)^2*diff(denseFLM$workGridY)[1],1,'sum')))
#  # densePredErr
  
  
#  # Errors are properly small
#  expect_lt(denseEstErr, 1) 
#  expect_gt(denseFLM$R2, 0.5) 
#  expect_lt(denseFLM$pv, 0.05) 
#  expect_equal(c(denseFLM$alpha), rep(alpha, M), tolerance = 0.2)
#  expect_equal(denseFLM$yHat, do.call(rbind, denseY$Ly), tolerance = 1)
#  expect_equal(denseFLM$yPred, do.call(rbind, denseYTest$Ly), tolerance = 1)
#  # expect_lt(densePredErr, 3) 

#  denseFLM$muY
#})

#test_that('Sparse, functional response case works', {
#  set.seed(1000)
  
#  library(MASS)
  
#  ### functional covariate
#  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
#  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
#  lambdaX <- c(1,0.7)
  
#  # training set
#  n <- 100
#  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
#  denseLt <- list(); denseLy <- list()
#  sparseLt <- list(); sparseLy <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:n) {
#    denseLt[[i]] <- t0
#    denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),3))
#    sparseLt[[i]] <- t0[ind]
#    sparseLy[[i]] <- denseLy[[i]][ind]
#  }
  
#  denseX <- list(Ly=denseLy,Lt=denseLt)
#  sparseX <- list(Ly=sparseLy,Lt=sparseLt)
  
#  denseX <- list(X=denseX)
#  sparseX <- list(X=sparseX)
  
#  # test set
#  N <- 500
  
#  XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
  
#  denseLtTest <- list(); denseLyTest <- list()
  
#  sparseLtTest <- list(); sparseLyTest <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:N) {
#    denseLtTest[[i]] <- t0
#    denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLtTest[[i]] <- t0[ind]
#    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#  }
  
#  denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#  sparseXTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
#  denseXTest <- list(X=denseXTest)
#  sparseXTest <- list(X=sparseXTest)
  
  
#  ### functional response
#  Beta <- matrix(c(1,0.5,-1,-1),nrow=2,ncol=2)
  
#  # training set
#  n <- 100
  
#  Eta <- Xi%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(n*2,0,0.5),nrow=n,ncol=2)
  
#  denseLt <- list()
#  denseLy <- list()
  
#  sparseLt <- list()
#  sparseLy <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:n) {
    
#    # comp 1
#    denseLt[[i]] <- t0
#    # denseLy[[i]] <- lambdaY[1]*Eta[i,1]*phi1(t0,1) + lambdaY[2]*Eta[i,2]*phi2(t0,1) + 
#    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
#    denseLy[[i]] <- Eta[i,1]*phi1(t0,1) + Eta[i,2]*phi2(t0,1) + 
#      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLt[[i]] <- t0[ind]
#    sparseLy[[i]] <- denseLy[[i]][ind]
#  }
  
#  denseY <- list(Ly=denseLy,Lt=denseLt)
#  sparseY <- list(Ly=sparseLy,Lt=sparseLt)
  
#  # test set
#  N <- 500
  
#  EtaTest <- XiTest%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(N*2,0,0.5),nrow=N,ncol=2) 
  
#  denseLtTest <- list()
#  denseLyTest <- list()
  
#  sparseLtTest <- list()
#  sparseLyTest <- list()
  
#  t0 <- seq(0,1,length.out=51)
#  for (i in 1:N) {
    
#    # comp 1
#    denseLtTest[[i]] <- t0
#    # denseLyTest[[i]] <- lambdaY[1]*EtaTest[i,1]*phi1(t0,1) + lambdaY[2]*EtaTest[i,2]*phi2(t0,1) + 
#    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
#    denseLyTest[[i]] <- EtaTest[i,1]*phi1(t0,1) + EtaTest[i,2]*phi2(t0,1) + 
#      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    
#    ind <- sort(sample(1:length(t0),5))
#    sparseLtTest[[i]] <- t0[ind]
#    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#  }
  
#  denseYTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#  sparseYTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
  
#  ## sparse
#  sparseFLM <- FLM1(Y=sparseY,X=sparseX,XTest=sparseXTest,optnsListX=list(FVEthreshold=0.95))
  
#  trueBetaList <- list()
#  trueBetaList[[1]] <- cbind(phi1(sparseFLM$workGridX[[1]],1),phi1(sparseFLM$workGridX[[1]],2))%*%Beta[1:2,]%*%
#    t(cbind(phi1(sparseFLM$workGridY,1),phi2(sparseFLM$workGridY,1)))
  
#  # coefficient function estimation error (L2-norm)
#  sparseEstErr <- sqrt(sum((sparseFLM$betaList[[1]] - trueBetaList[[1]])^2)*
#                         diff(sparseFLM$workGridX[[1]])[1]*
#                         diff(sparseFLM$workGridY)[1])
#  sparseEstErr
  
#  # par(mfrow=c(2,3))
#  # for (i in 1:6) {
#  #   i <- sample(1:n,1)
#  #   plot(sparseFLM$workGridY,sparseFLM$yHat[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (fitted)')
#  #   points(sparseY$Lt[[i]],sparseY$Ly[[i]],col=2)
#  # }
#  # 
#  # par(mfrow=c(2,3))
#  # for (i in 1:6) {
#  #   i <- sample(1:n,1)
#  #   plot(sparseFLM$workGridY,sparseFLM$yPred[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (pred)')
#  #   points(sparseYTest$Lt[[i]],sparseYTest$Ly[[i]],col=2)
#  # }
  
#  # # prediction error
#  # sparsePredErr <- sqrt(mean(apply((YTest - sparseFLM$yPred)^2*diff(sparseFLM$workGridY)[1],1,'sum')))
#  # sparsePredErr
  
  
  
#  # Errors are properly small
#  expect_lt(sparseEstErr, 3.61) 
#  # expect_lt(sparsePredErr, 3) 
#})


test_that('nRegGrid works', {
  set.seed(1000)
  
  library(MASS)
  
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  
  # training set
  n <- 100
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  denseLt <- list(); denseLy <- list()

  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLt[[i]] <- t0
    denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
  }
  
  denseX <- list(Ly=denseLy,Lt=denseLt)

  denseX <- list(X=denseX)

  # test set
  N <- 500
  
  XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
  
  denseLtTest <- list(); denseLyTest <- list()
  

  t0 <- seq(0,1,length.out=51)
  for (i in 1:N) {
    denseLtTest[[i]] <- t0
    denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
  }
  
  denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)

  denseXTest <- list(X=denseXTest)

  
  ### scalar response
  beta <- c(1, -1)
  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
  YTest <- c(XiTest%*%diag(lambdaX)%*%beta) + rnorm(N,0,0.5)
  
  ## dense
  denseFLM <- FLM1(Y=Y,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95,nRegGrid=100))
  
  expect_equal(length(denseFLM$workGridX[[1]]), 100)
})


test_that('A general test', {

  respType <- 'functional'
  predTypeCounts <- c(functional=1, scalar=1)
  funcType <- 'sparse' #TODO

for (ff in 0:1) {
for (ss in 0:2) {
  if (ff + ss == 0) next
  predTypeCounts <- c(functional = ff, scalar = ss)
  set.seed(1)

  n <- 50
  M <- 51
  K <- 2
  cBeta <- 1
  sigma <- 1
  lambdaY <- c(1, 1) # Must have length K
  basisType <- 'cos'

  domainX <- c(0, 1)
  domainY <- c(-0.5, 0.5)
  gridX <- seq(0, 1, length.out=M)
  gridY <- seq(domainY[1], domainY[2], length.out=M)

  BasisX <- CreateBasis(K, gridX, type=basisType)
  if (respType == 'scalar') {
    BasisY <- matrix(1)
    KY <- 1
  } else if (respType == 'functional') {
    BasisY <- BasisX
    KY <- K
  }

  muXFunc <- lapply(seq_len(predTypeCounts['functional']), 
                function(j) {function(x) x^(j -1) * 1})
  XFunc <- lapply(seq_len(predTypeCounts['functional']),
                  function(j) {
                    MakeSparseGP(n, muFun = muXFunc[[j]], sigma=sigma, basisType=basisType)
                  })
  # CreatePathPlot(inputData = XFunc[[2]])
  muXFuncGrid <- lapply(seq_along(XFunc), function(j) {
                    muXFunc[[j]](gridX)
                  })
  XTrue <- lapply(seq_along(XFunc), function(j) {
                    XFunc[[j]]$xi %*% t(BasisX) + 
                      matrix(muXFuncGrid[[j]], n, M, byrow=TRUE)
                  })
  muXScalar <- -as.numeric(seq_len(predTypeCounts['scalar']))
  XScalar <- lapply(seq_len(predTypeCounts['scalar']), 
                    function(j) {
                      muXScalar[j] + rnorm(n)
                    })
  BFunc <- matrix(rnorm(K * KY), K, KY)
  betaFunc <- lapply(seq_len(predTypeCounts['functional']),
                     function(j) {
                       cBeta * BasisX %*% BFunc %*% t(BasisY)
                     })
  betaScalar <- lapply(seq_len(predTypeCounts['scalar']), 
                       function(j) {
                         matrix(rnorm(K), nrow=1) %*% t(BasisY)
                       })
  beta <- c(betaFunc, betaScalar)
  betaEval <- list(do.call(rbind, betaFunc), 
                   do.call(rbind, betaScalar))
  betaEval <- betaEval[!vapply(betaEval, is.null, FALSE)]
  alpha <- gridY ^ 2 * 20

  if (respType == 'scalar') {
    alpha <- mean(alpha)
    MY <- 1
  } else {
    MY <- M
  }

  trueYGivenXZ <- matrix(alpha, n, MY, byrow=TRUE) + 
    Reduce(`+`, 
           mapply(function(X, beta) {
             X %*% beta * (gridX[2] - gridX[1])
           }, XTrue, betaFunc, SIMPLIFY=FALSE),
           init=matrix(0, n, MY)) +
    Reduce(`+`,
           mapply(function(Z, betaZ) {
             matrix(Z, ncol=1) %*% betaZ
           }, XScalar, betaScalar, SIMPLIFY=FALSE),
           init=matrix(0, n, MY))
  muY <- alpha + 
    Reduce(`+`, 
           mapply(function(X, beta) {
             X %*% beta * (gridX[2] - gridX[1])
           }, muXFuncGrid, betaFunc, SIMPLIFY=FALSE),
           init=rep(0, MY)) +
    Reduce(`+`,
           mapply(function(Z, betaZ) {
             matrix(Z, ncol=1) %*% betaZ
           }, muXScalar, betaScalar, SIMPLIFY=FALSE),
           init=rep(0, MY))
  muY <- c(muY)


  # matplot(t(trueYGivenXZ), type='l')

  # Add noise to Y
  if (respType == 'scalar') {
    YNoise <- rnorm(n, sd=sigma)
    Y <- trueYGivenXZ + YNoise
  } else if (respType == 'functional') {
    YNoise <- MakeSparseGP(n, sigma=sigma, lambda=lambdaY, basisType=basisType)
    YNoise$Lt <- lapply(YNoise$Lt, function(tt) tt * diff(domainY) + domainY[1])
    YNoise$Ly <- lapply(seq_len(n), function(i) {
                          tt <- YNoise$Lt[[i]]
                          yy <- YNoise$Ly[[i]] + 
                            ConvertSupport(gridY, tt, mu=trueYGivenXZ[i, ])
                          yy
                        })
    Y <- YNoise
  }

  res <- FLM1(Y, c(XFunc, XScalar), nPerm=1000)
  # b <- FLM(Y, XFunc, nPerm=1000)

  # TODO: add checks: optns passed in

  expect_equal(res$muY, muY, tolerance=0.5, scale=1)
  expect_equal(c(res$alpha), alpha, tolerance=1, scale=1)
  if (ff + ss < 3) {
    expect_equal(res$betaList, betaEval, tolerance=1, scale=1)
  } else {
    expect_equal(res$betaList, betaEval, tolerance=2, scale=1)
  }

  # # Plots
  # matplot(cbind(res$muY, muY), type='l')
  # matplot(cbind(c(res$alpha), alpha), type='l')
  # matplot(cbind(c(res$beta[[2]]), c(beta[[2]])), type='l')
  # matplot(t(res$betaList[[2]]), type='l', lty=1)
  # matplot(t(betaEval[[2]]), type='l', lty=2, add=TRUE)

  # image(res$betaList[[1]])
  # # image(res$betaList[[2]])
  # image(betaEval[[1]])
  # # image(betaFunc[[2]])

  # plot(c(res$alpha))
  # plot(alpha)
  # # matplot(t(res$betaList[[3]]), type='l', ylim=c(-5, 5))
  # matplot(t(do.call(rbind, betaScalar)), type='l', add=TRUE)

  # matplot(t(res$betaList[[1]]), type='l', ylim=c(-3, 3))
  # matplot(t(do.call(rbind, betaScalar)), type='l', add=TRUE)
  # matplot(t(res$yHat[1:5, ]), type='l')
  # matplot(t(trueYGivenXZ[1:5, ]), type='l', add=TRUE)

}
}
})
