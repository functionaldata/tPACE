# devtools::load_all()
set.seed(1)
n <- 100
pts <- seq(0, 1, by=0.05)
outPts <- seq(0, 1, by=0.1)
samp3 <- Wiener(n, pts) + rnorm(n * length(pts), sd=0.5)
samp3 <- Sparsify(samp3, pts, 5:10)
rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
brcov3 <- BinRawCov(rcov3)

gcv3b <- GCVLwls2DV2(pts, outPts, kern='epan', rcov=brcov3, t=samp3$Lt)

test_that('RotateLwls2dV2.R interface is correct', {
  expect_equal(Rrotatedmullwlsk(c(gcv3b$h, gcv3b$h) , 'epan', t(brcov3$tPairs), brcov3$meanVals, brcov3$count, rbind(outPts, outPts), npoly=1, bwCheck=FALSE), RotateLwls2DV2( gcv3b$h,  'epan', xin=brcov3$tPairs, yin=brcov3$meanVals, win=brcov3$count, xout=cbind(outPts, outPts)))
})
