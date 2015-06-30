devtools::load_all()
data(ethanol)

fitRef <- locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='epan')
tmp <- lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE)

# the 2D smoother does not match matlab version because it only uses elliptical window (rather than rectangular as in matlab PACE), and the 2D kernel implementations are different. The matlab epan kernel is not even the optimal spherical epan kernel, rather it is two epan kernels multiplied together.
test_that('The interface passes the arguments correctly', {
    expect_equal(fitted(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE, scale=TRUE), ethanol[, -1]), fitted(fitRef))
    expect_equal(fitted(lwls2d(0.5, kern='gauss', ethanol[, -1], ethanol[, 1], returnFit=TRUE, scale=TRUE), ethanol[, -1]), fitted(locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='gauss')))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], scale=TRUE), fitted(fitRef))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], xout=ethanol[1:5, 2:3], scale=TRUE), predict(fitRef, ethanol[1:5, 2:3]))
})

set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(n, pts, sparsify=length(pts))
mu3 <- rep(0, length(pts))
system.time(rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, 'Sparse', TRUE))
system.time(brcov3 <- BinRawCov(rcov3))
xout1 <- seq(0, 1, by=0.5)

test_that('Matlab version is similar to R version',
  expect_equal(as.numeric(lwls2d(0.3, kern='epan', rcov3$tPairs, rcov3$cxxn, xout1=xout1, xout2=xout1)), c(-0.0600,  0.0159,  0.0177,  0.0159,  0.4618,  0.5665,  0.0177,  0.5665,  1.1844), tolerance=0.1)
)

# Sparse case
set.seed(1)
n <- 400
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(n, pts, sparsify=5:10)
mu3 <- rep(0, length(pts))
system.time(rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, 'Sparse', TRUE))
system.time(brcov3 <- BinRawCov(rcov3))
xout1 <- seq(0, 1, by=0.1)


test_that('Estimated value should be close to the true value', {
  expect_equal(diag(lwls2d(0.3, kern='epan', rcov3$tPairs, rcov3$cxxn, xout1=xout1, xout2=xout1)), xout1, tolerance=0.15)
  expect_equal(diag(lwls2d(0.3, kern='epan', brcov3$tPairs, brcov3$meanVals, brcov3$count, xout1=xout1, xout2=xout1)), xout1, tolerance=0.15)
})

# Dense case
set.seed(1)
n <- 1000
pts <- seq(0, 1, by=0.01)
samp4 <- wiener(n, pts, sparsify=length(pts))
mu4 <- rep(0, length(pts))
system.time(rcov4 <- GetRawCov(samp4$yList, samp4$tList, pts, mu4, 'Dense', TRUE))
system.time(brcov4 <- BinRawCov(rcov4))
xout1 <- seq(0, 1, by=0.1)

test_that('Estimated value should be close to the true value', {
  expect_equal(diag(lwls2d(0.1, kern='epan', rcov4$tPairs, rcov4$cxxn, xout1=xout1, xout2=xout1)), xout1, tolerance=0.1)
  expect_equal(diag(lwls2d(0.1, kern='epan', brcov4$tPairs, brcov4$meanVals, brcov4$count, xout1=xout1, xout2=xout1)), xout1, tolerance=0.1)
})

xout2 <- c(0, 0.5, 1)
lwls2d(0.3, kern='rect', rcov3$tPairs, rcov3$cxxn, xout1=xout2, xout2=xout2)


# randWeights <- sample(1:3, n, replace=TRUE)
# library(R.matlab)
# writeMat('../../../data/data2Dsmoothers.mat', xin=rcov3$tPairs, yin=rcov3$cxxn, diagx=rcov3$diag[, 1], diagy=rcov3$diag[, 2], win=randWeights, xout1=xout1)
