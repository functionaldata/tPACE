devtools::load_all()
library(testthat)
# library(fda) # for generating fourier basis
# options(error=recover)
# trueLam <- 4 / ((2 * (1:50) - 1) * pi) ^ 2

set.seed(1)
n <- 500
p <- 500
pts <- seq(0, 1, length.out=p)
# samp <- matrix(rnorm(n), n, p)
samp <- wiener(n, pts)
samp <- samp + rnorm(n * p, sd=0.01) + matrix(pts, n, p, byrow=TRUE)
# matplot(t(samp[1:10, ]), type='l')
samp <- sparsify(samp, pts, p)
res <- FPCA(samp$yList, samp$tList, 
            list(dataType='Dense', error=FALSE, lean=TRUE))
# smRes <- FPCA(samp$yList, samp$tList, 
            # list(dataType='Dense', error=FALSE, lean=TRUE, 
                 # muCovEstMethod='smooth'))            
# fittedY <- fitted(res, derOptns=list(p=2))
# matplot(t(fittedY[1:10, ]), type='l', ylim=c(0, 2))


resDer <- FPCAder(res, list(p=1, GCV=TRUE))
# smDer <- FPCAder(smRes, list(p=1, GCV=TRUE))
plot(pts, resDer$mu, type='l', ylim=c(-1, 1))
plot(pts, resDer$muDer, type='l', ylim=c(-2, 2))
plot(pts, resDer$phi[, 1], type='l', ylim=c(-2, 2))
plot(pts, -resDer$phiDer[, 1], type='l')
plot(pts, resDer$phi[, 2], type='l', ylim=c(-2, 2))
plot(pts, -resDer$phiDer[, 2], type='l')
# plot(pts, smDer$mu, type='l', ylim=c(-1, 1))
# plot(pts, smDer$muDer, type='l')
# plot(pts, smDer$phi[, 1], type='l', ylim=c(-2, 2))
# plot(pts, smDer$phiDer[, 1], type='l')
test_that('first derivative is fine', {
  expect_equal(mean(resDer$muDer), 1, tolerance=0.1)
})


# test_that('FPCAder returns correct derivatives of mean and eigenfunctions for both dense and sparse cases',{
  # set.seed(1)
  # n <- 100
  # pts <- seq(0, 1, by=0.01)
  # #mu <- rep(0, length(pts))
  # mu <- 0:(length(pts)-1) / 50
  # samp1 <- t(mu + fourier(pts, nbasis=length(trueLam))[,1:50] %*% matrix(rnorm(n*length(trueLam),mean=0,sd=sqrt(rep(trueLam,n))),nrow=length(trueLam))+
    # rnorm(n * length(pts), mean = 0, sd = 0.1))
  # samp2 <- sparsify(samp1, pts, 10)
  # samp11 <- list(tList = c(), yList = c())
  # for(i in 1:n){samp11$tList[[i]] <- pts; samp11$yList[[i]] <- samp1[i,]}
  # samp1 <- samp11
  # p1 <- SetOptions(samp1$yList, samp1$tList, list(dataType='Dense', error=TRUE, kernel='epan', verbose=TRUE))
  # p2 <- SetOptions(samp2$yList, samp2$tList, list(dataType='Sparse', error=TRUE, kernel='epan', verbose=TRUE))
  # set.seed(1)
  # FPCAres1 <- FPCA(samp1$yList, samp1$tList, p1)
  # Derres1 <-  FPCAder(FPCAres1)
  # set.seed(1)
  # FPCAres2 <- FPCA(samp2$yList, samp2$tList, p2)
  # Derres2 <- FPCAder(FPCAres2)
  
  # # check equality of derivatives of mean and eigenfunctions
  # #expect_equal(Derres1$mu, rep(0.02, length(Derres1$mu)))
  # expect_equal(Derres1$mu, getDerivative(y=FPCAres1$mu, t=FPCAres1$workGrid))
  # #expect_equal(Derres2$mu, rep(0.02, length(Derres2$mu)))
  # expect_equal(Derres2$mu, getDerivative(y=FPCAres2$mu, t=FPCAres2$obsGrid))  
  # #expect_equal(Derres1$phi, fourier(FPCAres1$workGrid,nbasis=length(trueLam),nderiv=1)[,1:ncol(Derres1$phi)])
  # expect_equal(Derres1$phi, apply(FPCAres1$phi,2,getDerivative,t=FPCAres1$workGrid))
  # #expect_equal(Derres2$phi, fourier(FPCAres2$workGrid,nbasis=length(trueLam),nderiv=1)[,1:ncol(Derres2$phi)])
  # expect_equal(Derres2$phi, apply(FPCAres2$phi,2,getDerivative,t=FPCAres2$workGrid))
# })
