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

test_that('noisy dense case for DPC', {
  bw <- 0.1
  set.seed(1)
  n <- 100 
  plotInd <- 1:8
  pts <- seq(0, 1, by=0.01)
  sparsity <- length(pts)
  mu <- seq(0, 1, length.out=length(pts))
  phi <-  CreateBasis(K=3, type='fourier',  pts= pts)
  lambdaTrue <- c(1, 0.8, 0.5)^2
  sigma2 <- 1
  nRegGrid <- 51
  regGrid <- seq(0, 1, length.out=nRegGrid)

  samp1 <- t(mu + phi %*% matrix(rnorm(n*3, mean=0, sd=sqrt(lambdaTrue)), nrow=3))
  samp1p <- t(apply(samp1, 1, diff) / (pts[2] - pts[1]))
  samp1p <- t(ConvertSupport(seq(min(pts), max(pts), length.out=length(pts) - 1), regGrid, phi=t(samp1p)))
  samp2 <- Sparsify(samp1 + rnorm(length(samp1), sd=sqrt(sigma2)), pts, sparsity)
  fpcaObj <- FPCA(samp2$Lt, Ly= samp2$Ly, list(nRegGrid=nRegGrid, methodMuCovEst='smooth', userBwCov=bw, userBwMu=bw)) 

  derOptns1 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=TRUE)
  fpcaObjDer1 <- FPCAder(fpcaObj, derOptns1)
  derOptns2 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=FALSE)
  fpcaObjDer2 <- FPCAder(fpcaObj, derOptns2)
  k1 <- 2
  est1 <- fitted(fpcaObjDer1, k1)
  k2 <- 2
  est2 <- fitted(fpcaObjDer2, k2)
  expect_equal(est1, samp1p, tolerance=2, scale=1)
# 1D smoother for cov10 is better
  expect_true(mean(abs(est1 - samp1p)) < mean(abs(est2 - samp1p))) 

# # display.brewer.all()
# color <- brewer.pal(8, 'Dark2')
# ylim <- c(min(samp1p[, plotInd]), max(samp1p[, plotInd]))
# matplot(regGrid, t(samp1p[plotInd, ]), type='l', ylim=ylim, col=color, lwd=2)
# CreatePathPlot(fpcaObjDer1, plotInd, k1, derOptns=derOptns1, ylim=ylim, col=color, lwd=2)
# CreatePathPlot(fpcaObjDer2, plotInd, k2, derOptns=derOptns2, ylim=ylim, col=color, lwd=2)
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

# set.seed(1)
# n <- 100 
# plotInd <- 1:8
# pts <- seq(0, 1, by=0.01)
# sparsity <- length(pts)
# # mu <- 0:(length(pts)-1) / 50
# mu <- rep(0, length(pts))
# phi <-  CreateBasis(K=3, type='fourier',  pts= pts)
# lambdaTrue <- c(1, 0.8, 0.5)^2
# covTrue <- phi %*% diag(lambdaTrue) %*% t(phi)
# # covTrue1 <- sapply(pts, function(t) sapply(pts, function(s) 1 + 0.5 * sin(2 * pi * s) * sin(2 * pi * t)))
# # expect_equal(covTrue, covTrue1)
# cov10True <- sapply(pts, function(t) sapply(pts, function(s) 4 * pi * lambdaTrue[2] * cos(2 * pi * s) * sin(2 * pi * t) - 4 * pi * lambdaTrue[3] * sin(2 * pi * s) * cos(2 * pi * t)))
# cov11True <- sapply(pts, function(t) sapply(pts, function(s) 8 * pi^2 * lambdaTrue[2] * cos(2 * pi * s) * cos(2 * pi * t) + 8 * pi^2 * lambdaTrue[3] * sin(2 * pi * s) * cos(2 * pi * t)))
# # cov10True1 <- apply(covTrue, 2, diff) / (pts[2] - pts[1])
# # cov11True1 <- apply(apply(covTrue, 2, diff), 1, diff) / (pts[2] - pts[1])^2
# rgl::persp3d(pts, pts, covTrue, xlab='s', ylab='t')
# rgl::persp3d(pts, pts, cov10True, xlab='s', ylab='t')
# rgl::persp3d(pts, pts, cov11True, xlab='s', ylab='t')
# # rgl::persp3d(pts[-1], pts, cov10True1, xlab='s', ylab='t')
# # rgl::persp3d(pts[-1], pts[-1], cov11True1)

# samp1 <- t(mu + phi %*% matrix(rnorm(n*3, mean=0, sd=sqrt(lambdaTrue)), nrow=3))
# samp1p <- apply(samp1, 1, diff) / (pts[2] - pts[1])
# matplot(pts, t(samp1[plotInd, ]), type='l')
# samp2 <- Sparsify(samp1, pts, sparsity)
# fpcaObj <- FPCA(samp2$Lt, Ly= samp2$Ly, list(nRegGrid=51))
# derOptns <- list(bw=0.1, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=TRUE)
# fpcaObjDer <- FPCAder(fpcaObj, derOptns)

# # display.brewer.all()
# color <- brewer.pal(8, 'Dark2')
# ylim <- c(min(samp1p[, plotInd]), max(samp1p[, plotInd]))
# matplot(pts[-1], samp1p[, plotInd], type='l', ylim=ylim, col=color, lwd=2)
# CreatePathPlot(fpcaObjDer, plotInd, min(length(lambdaTrue), length(fpcaObjDer$lambdaDer)), derOptns=derOptns, ylim=ylim, col=color, lwd=2)

test_that('noisy sparse case for DPC', {
  set.seed(1)
  bw <- 0.2
  n <- 200 
  plotInd <- 1:8
  pts <- seq(0, 1, by=0.01)
  sparsity <- 2:9
  mu <- seq(0, 1, length.out=length(pts))
  phi <-  CreateBasis(K=3, type='fourier',  pts= pts)
  lambdaTrue <- c(1, 0.8, 0.5)^2
  sigma2 <- 1
  nRegGrid <- 51
  regGrid <- seq(0, 1, length.out=nRegGrid)

  samp1 <- t(mu + phi %*% matrix(rnorm(n*3, mean=0, sd=sqrt(lambdaTrue)), nrow=3))
  samp1p <- t(apply(samp1, 1, diff) / (pts[2] - pts[1]))
  samp1p <- t(ConvertSupport(seq(min(pts), max(pts), length.out=length(pts) - 1), regGrid, phi=t(samp1p)))
  samp2 <- Sparsify(samp1 + rnorm(length(samp1), sd=sqrt(sigma2)), pts, sparsity)
  fpcaObj <- FPCA(samp2$Lt, Ly= samp2$Ly, list(nRegGrid=nRegGrid, methodMuCovEst='smooth', userBwCov=bw, userBwMu=bw)) 

  derOptns1 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=TRUE)
  fpcaObjDer1 <- FPCAder(fpcaObj, derOptns1)
  derOptns2 <- list(bw=bw, kernelType='epan', p=1, method='DPC', useTrue=FALSE, G10_1D=FALSE)
  fpcaObjDer2 <- FPCAder(fpcaObj, derOptns2)
  k1 <- 2
  est1 <- fitted(fpcaObjDer1, k1)
  k2 <- 2
  est2 <- fitted(fpcaObjDer2, k2)
  CreatePathPlot(fpcaObjDer1, plotInd, ylim=c(-15, 20))
  CreatePathPlot(fpcaObjDer2, plotInd, ylim=c(-15, 20))
  matplot(t(samp1p[plotInd, ]), type='l', ylim=c(-15, 20))

  expect_equal(est1, samp1p, tolerance=5, scale=1)
  expect_equal(est2, samp1p, tolerance=5, scale=1)
})
