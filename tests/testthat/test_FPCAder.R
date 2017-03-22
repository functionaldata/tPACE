# devtools::load_all()
library(testthat)
# library(RColorBrewer)

test_that('FPCAder correct derivatives of mean for dense case', {
  set.seed(1)
  n <- 500
  p <- 300
  pts <- seq(0, 1, length.out=p)
  samp <- Wiener(n, pts)
  samp <- samp + rnorm(n * p, sd=0.01) + matrix(pts, n, p, byrow=TRUE)
  spSamp <- Sparsify(samp, pts, p) # This just puts the sample in list.
  fpcaObj <- FPCA(spSamp$Ly, spSamp$Lt, list(dataType='Dense' ))
  fpcaObjDer <- FPCAder(fpcaObj, list(p=1))
    
  expect_equal(median(fpcaObjDer$muDer), 1, tolerance=0.1, scale = 1)
})




test_that('FPCAder correct derivatives of mean for sparse case',{
  set.seed(1)
  n <- 333 
  pts <- seq(0, 1, by=0.01)
  mu <- 0:(length(pts)-1) / 50
  phi =  CreateBasis(K=4, type='fourier',  pts= pts)
  samp1 <- t(mu + phi %*% matrix( rnorm(n*4, mean=0, sd=c(1,.4, 0.01, 0.001)), nrow=4))
  samp2 <- Sparsify(samp1, pts, 8)
  fpcaObj = FPCA(samp2$Lt, Ly= samp2$Ly)
  fpcaObjDer <- FPCAder(fpcaObj, list(p=1))
  expect_equal(median(fpcaObjDer$muDer), 2, tolerance=0.1, scale = 1) 
})


test_that('noisy dense case for DPC, FPC, and FPC1', {
  bw <- 0.2
  kern <- 'epan'
  set.seed(1)
  n <- 100 
  plotInd <- 1:8
  pts <- seq(0, 1, by=0.01)
  M <- length(pts)
  sparsity <- length(pts)
  mu <- seq(0, 1, length.out=length(pts))
  phi <-  CreateBasis(K=3, type='legendre01', pts=pts)
  basisDerMat <- apply(phi, 2, function(x) ConvertSupport(seq(0, 1, length.out=M - 1), pts, diff(x) * (M - 1)))
  lambdaTrue <- c(1, 0.8, 0.5)^2
  sigma2 <- 1

  samp2 <- MakeGPFunctionalData(n, M, pts, K=length(lambdaTrue), lambda=lambdaTrue, sigma=sqrt(sigma2), basisType='legendre01')
  samp2 <- c(samp2, MakeFPCAInputs(tVec=pts, yVec=samp2$Yn))
  trueDer <- matrix(1, n, M, byrow=TRUE) + tcrossprod(samp2$xi, basisDerMat)
  samp2 <- c(samp2, list(trueDer=trueDer))
  fpcaObj <- FPCA(samp2$Ly, samp2$Lt, list(nRegGrid=M, methodMuCovEst='smooth', userBwCov=bw, userBwMu=bw, kernel=kern)) 

  FPCoptn1 <- list(bw=bw, kernelType=kern, p=1, method='FPC')
  FPCoptn2 <- list(bw=bw, kernelType=kern, p=1, method='FPC1')
  DPCoptn1 <- list(bw=bw, kernelType=kern, p=1, method='DPC', G10_1D=TRUE)
  DPCoptn2 <- list(bw=bw, kernelType=kern, p=1, method='DPC', G10_1D=FALSE)
  FPC <- FPCAder(fpcaObj, FPCoptn1)
  FPC1 <- FPCAder(fpcaObj, FPCoptn2)
  fpcaObjDer1 <- FPCAder(fpcaObj, DPCoptn1)
  fpcaObjDer2 <- FPCAder(fpcaObj, DPCoptn2)
  k1 <- 3
  estFPC <- fitted(FPC, k1)
  estFPC1 <- fitted(FPC1, k1)
  estDPC <- fitted(fpcaObjDer1, k1)
  estDPC1 <- fitted(fpcaObjDer2, k1)

  matplot(t(samp2$Y), type='l')
  matplot(t(trueDer), type='l')
  matplot(t(estFPC1), type='l')
  matplot(t(estFPC), type='l')
  matplot(t(estDPC), type='l')
  matplot(t(estDPC1), type='l')

  expect_equal(estFPC, trueDer, tolerance=1, scale=1)
  expect_equal(estFPC1, trueDer, tolerance=1, scale=1)
  expect_equal(estDPC, trueDer, tolerance=1, scale=1)

  # 1D smoother for cov10 is better
  expect_true(mean(abs(estDPC - trueDer)) < mean(abs(estDPC1 - trueDer))) 
})

# test_that('noiseless dense case for DPC', {
  # bw <- 0.1
  # set.seed(1)
  # n <- 100 
  # plotInd <- 1:8
  # pts <- seq(0, 1, by=0.01)
  # sparsity <- length(pts)
  # mu <- seq(0, 1, length.out=length(pts))
  # phi <-  CreateBasis(K=3, type='fourier',  pts= pts)
  # lambdaTrue <- c(1, 0.8, 0.5)^2
  # sigma2 <- 0
  # nRegGrid <- 51
  # regGrid <- seq(0, 1, length.out=nRegGrid)

  # samp1 <- t(mu + phi %*% matrix(rnorm(n*3, mean=0, sd=sqrt(lambdaTrue)), nrow=3))
  # samp1p <- t(apply(samp1, 1, diff) / (pts[2] - pts[1]))
  # samp1p <- t(ConvertSupport(seq(min(pts), max(pts), length.out=length(pts) - 1), regGrid, phi=t(samp1p)))
  # samp2 <- Sparsify(samp1 + rnorm(length(samp1), sd=sqrt(sigma2)), pts, sparsity)
  # # Cross-sectional is used!
  # fpcaObj <- FPCA(samp2$Lt, Ly= samp2$Ly, list(nRegGrid=nRegGrid, methodMuCovEst='cross-sectional'))

  # derOptns1 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=TRUE)
  # fpcaObjDer1 <- FPCAder(fpcaObj, derOptns1)
  # derOptns2 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=FALSE)
  # fpcaObjDer2 <- FPCAder(fpcaObj, derOptns2)
  # k1 <- 2
  # est1 <- fitted(fpcaObjDer1, k1, derOptns1)
  # k2 <- 2
  # est2 <- fitted(fpcaObjDer2, k2, derOptns2)
  # expect_equal(est1, samp1p, tolerance=1.5, scale=1)
  # # 1D smoother for cov10 is better
  # expect_true(mean(abs(est1 - samp1p)) < mean(abs(est2 - samp1p))) 

  # # # display.brewer.all()
  # # color <- brewer.pal(8, 'Dark2')
  # # ylim <- c(min(samp1p[, plotInd]), max(samp1p[, plotInd]))
  # # matplot(regGrid, t(samp1p[plotInd, ]), type='l', ylim=ylim, col=color, lwd=2)
  # # CreatePathPlot(fpcaObjDer1, plotInd, 2, derOptns=derOptns1, ylim=ylim, col=color, lwd=2)
  # # CreatePathPlot(fpcaObjDer2, plotInd, 2, derOptns=derOptns2, ylim=ylim, col=color, lwd=2)
# })


test_that('noisy sparse case for DPC, FPC, and FPC1', {
  bw <- 0.4
  kern <- 'epan'
  set.seed(1)
  n <- 100 
  pts <- seq(0, 1, by=0.01)
  M <- length(pts)
  sparsity <- 2:9
  mu <- seq(0, 1, length.out=length(pts))
  phi <-  CreateBasis(K=3, type='legendre01', pts=pts)
  basisDerMat <- apply(phi, 2, function(x) ConvertSupport(seq(0, 1, length.out=M - 1), pts, diff(x) * (M - 1)))
  lambdaTrue <- c(1, 0.8, 0.5)^2
  sigma2 <- 0.01

  samp2 <- MakeSparseGP(n, runif, sparsity, identity, K=length(lambdaTrue), lambda=lambdaTrue, sigma=sqrt(sigma2), basisType='legendre01')
  trueDer <- matrix(1, n, M, byrow=TRUE) + tcrossprod(samp2$xi, basisDerMat)
  samp2 <- c(samp2, list(trueDer=trueDer))
  fpcaObj <- FPCA(samp2$Ly, samp2$Lt, list(nRegGrid=M, methodMuCovEst='smooth', userBwCov=bw, userBwMu=bw, kernel=kern)) 

  FPCoptn1 <- list(bw=bw, kernelType=kern, p=1, method='FPC')
  FPCoptn2 <- list(bw=bw, kernelType=kern, p=1, method='FPC1')
  DPCoptn1 <- list(bw=bw, kernelType=kern, p=1, method='DPC', G10_1D=TRUE)
  DPCoptn2 <- list(bw=bw, kernelType=kern, p=1, method='DPC', G10_1D=FALSE)
  FPC <- FPCAder(fpcaObj, FPCoptn1)
  FPC1 <- FPCAder(fpcaObj, FPCoptn2)
  fpcaObjDer1 <- FPCAder(fpcaObj, DPCoptn1)
  fpcaObjDer2 <- FPCAder(fpcaObj, DPCoptn2)
  k1 <- 2
  estFPC <- fitted(FPC, k1)
  estFPC1 <- fitted(FPC1, k1)
  estDPC <- fitted(fpcaObjDer1, k1)
  estDPC1 <- fitted(fpcaObjDer2, k1)

  matplot(t(trueDer), type='l')
  matplot(t(estFPC), type='l')
  matplot(t(estFPC1), type='l')
  matplot(t(estDPC), type='l')
  matplot(t(estDPC1), type='l')

  expect_equal(estFPC, trueDer, tolerance=5, scale=1)
  expect_equal(estFPC1, trueDer, tolerance=5, scale=1)
  expect_equal(estDPC, trueDer, tolerance=2, scale=1)

  # 1D smoother for cov10 is better
  expect_true(mean(abs(estDPC - trueDer)) < mean(abs(estDPC1 - trueDer))) 
})
