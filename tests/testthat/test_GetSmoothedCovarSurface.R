devtools::load_all()
library(testthat)

# GMeanAndGCV

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
samp3 <- wiener(50, pts, sparsify=length(pts))
mu3 <- rep(0, length(pts))

# without error
p0 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(dataType='Sparse', error=FALSE, kernel='epan'))
noErr <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p0, useBins=FALSE)

# with error
p1 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(dataType='Sparse', error=TRUE, kernel='epan'))
Err <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p1, useBins=FALSE)

# unit tests: test the interface.
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, p1$dataType, p1$error)
system.time(tmp <- gcvlwls2dV2(pts, regGrid=regGrid, kern='epan', rcov=rcov3, t=samp3$tList))
gcvBW3 <- sqrt(tmp$h * tmp$minBW)
sigma23 <- pc_covE(pts, regGrid, gcvBW3, kernel='epan', rcov=rcov3, rotationCut=c(0.25, 0.75))$sigma2
test_that('Smooth Cov Surface interface is right', {
  expect_equal(rcov3, Err$rawCov)
  expect_equal(gcvBW3, Err$bwCov)
  expect_equal(sigma23, Err$sigma2)
})

# GCV
p2 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(bwuserCovGcv='GCV', dataType='Sparse', error=FALSE, kernel='epan'))
tmp2 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p2, useBins=FALSE)
sum((diag(tmp2$smoothCov) - seq(0, 1, by=0.1))^2)

# CV
p3 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(bwuserCovGcv='CV', dataType='Sparse', error=FALSE, kernel='epan'))
system.time(tmp3 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p3, useBins=FALSE))
sum((diag(tmp3$smoothCov) - seq(0, 1, by=0.1))^2)

# truncation.
pTrunc <- SetOptions(samp3$yList, samp3$tList, CreateOptions(dataType='Sparse', error=FALSE, kernel='epan', outPercent=c(0.01, 0.99)))
noErrTrunc <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, pTrunc, useBins=FALSE)
test_that('Cov Surface truncation works', {
  expect_equal(noErr$smoothedCov[2:10, 2:10], noErrTrunc$smoothedCov)
  expect_equal(noErr$rawCov, noErrTrunc$rawCov)
  expect_equal(noErr$bwCov, noErrTrunc$bwCov)
  expect_equal(noErr$sigma2, noErrTrunc$sigma2)
})
