library(MASS)
library(testthat)
#devtools::load_all()

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
  denseFLM <- FLM(Y=Y,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95))
  
  trueBetaList <- list()
  trueBetaList[[1]] <- cbind(phi1(denseFLM$workGridX[[1]],1),phi1(denseFLM$workGridX[[1]],2))%*%beta
  
  # coefficient function estimation error (L2-norm)
  par(mfrow=c(1,1))
  plot(denseFLM$workGridX[[1]],denseFLM$betaList[[1]],type='l',xlab='t',ylab=paste('beta',1,sep=''))
  points(denseFLM$workGridX[[1]],trueBetaList[[1]],type='l',col=2)
  
  denseEstErr <- sqrt(trapzRcpp(denseFLM$workGridX[[1]],(denseFLM$betaList[[1]] - trueBetaList[[1]])^2))
  denseEstErr
  
  par(mfrow=c(1,2))
  plot(denseFLM$yHat,Y,xlab='fitted Y', ylab='observed Y')
  abline(coef=c(0,1),col=8)
  
  plot(denseFLM$yPred,YTest,xlab='predicted Y', ylab='observed Y')
  abline(coef=c(0,1),col=8)
  
  # prediction error
  densePredErr <- sqrt(mean((YTest - denseFLM$yPred)^2))
  densePredErr
  
  # Errors are properly small
  expect_lt(denseEstErr, 0.1) 
  expect_lt(densePredErr, 1) 
})



test_that('Sparse, scalar response case works', {
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
  
  sparseFLM <- FLM(Y=Y,X=sparseX,XTest=sparseXTest,optnsListX=list(FVEthreshold=0.95))
  
  trueBetaList <- list()
  trueBetaList[[1]] <- cbind(phi1(sparseFLM$workGridX[[1]],1),phi1(sparseFLM$workGridX[[1]],2))%*%beta
  
  # coefficient function estimation error (L2-norm)
  par(mfrow=c(1,1))
  plot(sparseFLM$workGridX[[1]],sparseFLM$betaList[[1]],type='l',xlab='t',ylab=paste('beta',1,sep=''))
  points(sparseFLM$workGridX[[1]],trueBetaList[[1]],type='l',col=2)
  
  sparseEstErr <- sqrt(trapzRcpp(sparseFLM$workGridX[[1]],(sparseFLM$betaList[[1]] - trueBetaList[[1]])^2))
  sparseEstErr
  
  par(mfrow=c(1,2))
  plot(sparseFLM$yHat,Y,xlab='fitted Y', ylab='observed Y')
  abline(coef=c(0,1),col=8)
  
  plot(sparseFLM$yPred,YTest,xlab='fitted Y', ylab='observed Y')
  abline(coef=c(0,1),col=8)
  
  # prediction error
  sparsePredErr <- sqrt(mean((YTest - sparseFLM$yPred)^2))
  sparsePredErr
  
  
  # Errors are properly small
  expect_lt(sparseEstErr, 3.1) 
  expect_lt(sparsePredErr, 3) 
})





test_that('Dense, functional response case works', {
  
  M <- 51
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
  
  t0 <- seq(0,1,length.out=M)
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
  
  t0 <- seq(0,1,length.out=M)
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
  
  
  ### functional response
  alpha <- 1
  Beta <- matrix(c(1,0.5,-1,-1),nrow=2,ncol=2)

  # training set
  n <- 100
  
  Eta <- Xi%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(n*2,0,0.5),nrow=n,ncol=2)
  
  denseLt <- list()
  denseLy <- list()
  
  sparseLt <- list()
  sparseLy <- list()
  
  t0 <- seq(0,1,length.out=M)
  for (i in 1:n) {
    
    # comp 1
    denseLt[[i]] <- t0
    # denseLy[[i]] <- lambdaY[1]*Eta[i,1]*phi1(t0,1) + lambdaY[2]*Eta[i,2]*phi2(t0,1) + 
    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    denseLy[[i]] <- Eta[i,1]*phi1(t0,1) + Eta[i,2]*phi2(t0,1) + 
      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) + alpha
    
    ind <- sort(sample(1:length(t0),5))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[[i]][ind]
  }
  
  denseY <- list(Ly=denseLy,Lt=denseLt)
  sparseY <- list(Ly=sparseLy,Lt=sparseLt)
  
  # test set
  N <- 500
  
  EtaTest <- XiTest%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(N*2,0,0.5),nrow=N,ncol=2) 
  
  denseLtTest <- list()
  denseLyTest <- list()
  
  sparseLtTest <- list()
  sparseLyTest <- list()
  
  t0 <- seq(0,1,length.out=M)
  for (i in 1:N) {
    
    # comp 1
    denseLtTest[[i]] <- t0
    # denseLyTest[[i]] <- lambdaY[1]*EtaTest[i,1]*phi1(t0,1) + lambdaY[2]*EtaTest[i,2]*phi2(t0,1) + 
    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    denseLyTest[[i]] <- EtaTest[i,1]*phi1(t0,1) + EtaTest[i,2]*phi2(t0,1) + 
      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) + alpha
    
    ind <- sort(sample(1:length(t0),5))
    sparseLtTest[[i]] <- t0[ind]
    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
  }
  
  denseYTest <- list(Ly=denseLyTest,Lt=denseLtTest)
  sparseYTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
  
  ## dense
  # print(system.time({
  denseFLM <- FLM(Y=denseY,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95), nPerm=1000)
  # }))
  
  trueBetaList <- list()
  trueBetaList[[1]] <- cbind(phi1(denseFLM$workGridX[[1]],1),phi1(denseFLM$workGridX[[1]],2))%*%Beta[1:2,]%*%
    t(cbind(phi1(denseFLM$workGridY,1),phi2(denseFLM$workGridY,1)))
  
  # coefficient function estimation error (L2-norm)
  denseEstErr <- sqrt(sum((denseFLM$betaList[[1]] - trueBetaList[[1]])^2)*
                        diff(denseFLM$workGridX[[1]])[1]*
                        diff(denseFLM$workGridY)[1])
  denseEstErr
  
  # par(mfrow=c(2,3))
  # for (i in 1:6) {
  #   i <- sample(1:n,1)
  #   plot(denseFLM$workGridY,denseFLM$yHat[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (hat)')
  #   points(denseY$Lt[[i]],denseY$Ly[[i]],type='l',col=2)
  # }
  # 
  # par(mfrow=c(2,3))
  # for (i in 1:6) {
  #   i <- sample(1:n,1)
  #   plot(denseFLM$workGridY,denseFLM$yPred[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (pred)')
  #   points(denseYTest$Lt[[i]],denseYTest$Ly[[i]],type='l',col=2)
  # }
  
  # # prediction error
  # densePredErr <- sqrt(mean(apply((YTest - denseFLM$yPred)^2*diff(denseFLM$workGridY)[1],1,'sum')))
  # densePredErr
  
  
  # Errors are properly small
  expect_lt(denseEstErr, 1) 
  expect_gt(denseFLM$R2, 0.5) 
  expect_lt(denseFLM$pv, 0.05) 
  expect_equal(denseFLM$alpha, rep(alpha, M), tolerance = 0.2)
  expect_equal(denseFLM$yHat, do.call(rbind, denseY$Ly), tolerance = 1)
  expect_equal(denseFLM$yPred, do.call(rbind, denseYTest$Ly), tolerance = 1)
  # expect_lt(densePredErr, 3) 
})

test_that('Sparse, functional response case works', {
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
  
  
  ### functional response
  Beta <- matrix(c(1,0.5,-1,-1),nrow=2,ncol=2)
  
  # training set
  n <- 100
  
  Eta <- Xi%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(n*2,0,0.5),nrow=n,ncol=2)
  
  denseLt <- list()
  denseLy <- list()
  
  sparseLt <- list()
  sparseLy <- list()
  
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    
    # comp 1
    denseLt[[i]] <- t0
    # denseLy[[i]] <- lambdaY[1]*Eta[i,1]*phi1(t0,1) + lambdaY[2]*Eta[i,2]*phi2(t0,1) + 
    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    denseLy[[i]] <- Eta[i,1]*phi1(t0,1) + Eta[i,2]*phi2(t0,1) + 
      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    
    ind <- sort(sample(1:length(t0),5))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[[i]][ind]
  }
  
  denseY <- list(Ly=denseLy,Lt=denseLt)
  sparseY <- list(Ly=sparseLy,Lt=sparseLt)
  
  # test set
  N <- 500
  
  EtaTest <- XiTest%*%diag(lambdaX)%*%Beta #+ matrix(rnorm(N*2,0,0.5),nrow=N,ncol=2) 
  
  denseLtTest <- list()
  denseLyTest <- list()
  
  sparseLtTest <- list()
  sparseLyTest <- list()
  
  t0 <- seq(0,1,length.out=51)
  for (i in 1:N) {
    
    # comp 1
    denseLtTest[[i]] <- t0
    # denseLyTest[[i]] <- lambdaY[1]*EtaTest[i,1]*phi1(t0,1) + lambdaY[2]*EtaTest[i,2]*phi2(t0,1) + 
    #   rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    denseLyTest[[i]] <- EtaTest[i,1]*phi1(t0,1) + EtaTest[i,2]*phi2(t0,1) + 
      rnorm(1,0,0.25)*phi1(t0,2) + rnorm(1,0,0.25)*phi2(t0,2) 
    
    ind <- sort(sample(1:length(t0),5))
    sparseLtTest[[i]] <- t0[ind]
    sparseLyTest[[i]] <- denseLyTest[[i]][ind]
  }
  
  denseYTest <- list(Ly=denseLyTest,Lt=denseLtTest)
  sparseYTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
  
  
  ## sparse
  sparseFLM <- FLM(Y=sparseY,X=sparseX,XTest=sparseXTest,optnsListX=list(FVEthreshold=0.95))
  
  trueBetaList <- list()
  trueBetaList[[1]] <- cbind(phi1(sparseFLM$workGridX[[1]],1),phi1(sparseFLM$workGridX[[1]],2))%*%Beta[1:2,]%*%
    t(cbind(phi1(sparseFLM$workGridY,1),phi2(sparseFLM$workGridY,1)))
  
  # coefficient function estimation error (L2-norm)
  sparseEstErr <- sqrt(sum((sparseFLM$betaList[[1]] - trueBetaList[[1]])^2)*
                         diff(sparseFLM$workGridX[[1]])[1]*
                         diff(sparseFLM$workGridY)[1])
  sparseEstErr
  
  # par(mfrow=c(2,3))
  # for (i in 1:6) {
  #   i <- sample(1:n,1)
  #   plot(sparseFLM$workGridY,sparseFLM$yHat[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (fitted)')
  #   points(sparseY$Lt[[i]],sparseY$Ly[[i]],col=2)
  # }
  # 
  # par(mfrow=c(2,3))
  # for (i in 1:6) {
  #   i <- sample(1:n,1)
  #   plot(sparseFLM$workGridY,sparseFLM$yPred[i,],type='l',ylim=c(-10,10),xlab='t',ylab='Y (pred)')
  #   points(sparseYTest$Lt[[i]],sparseYTest$Ly[[i]],col=2)
  # }
  
  # # prediction error
  # sparsePredErr <- sqrt(mean(apply((YTest - sparseFLM$yPred)^2*diff(sparseFLM$workGridY)[1],1,'sum')))
  # sparsePredErr
  
  
  
  # Errors are properly small
  expect_lt(sparseEstErr, 3.61) 
  # expect_lt(sparsePredErr, 3) 
})


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
  denseFLM <- FLM(Y=Y,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95,nRegGrid=100))
  
  expect_equal(length(denseFLM$workGridX[[1]]), 100)
})


