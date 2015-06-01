devtools::load_all()
options(error=recover)
library(testthat)

load('../data/dataForGetRawCov.RData')
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output

# getGCVscores
test_that('getGCVscores interface works', {
    expect_equal(
        getGCVscores(2, kern='epan', xin=rcov$tpairn, yin=rcov$cxxn, win=as.numeric(rcov$win)), 
        unname(gcv(rcov$cxxn ~ lp(rcov$tpairn[, 1], rcov$tpairn[, 2], deg=1, h=2), kern='epan')['gcv'])
    )
})


# use wiener process
set.seed(1)
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(100, pts, sparsify=20)
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
system.time(h3 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=rcov3))
fit3 <- lwls2d(h3$h, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, returnFit=TRUE)
plot(fit3)

# Binned rcov = rcov
brcov3 <- BinRawCov(rcov3)
system.time(h3 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=rcov3))
system.time(h3bin <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=brcov3))

# Error=TRUE
samp4 <- samp3
samp4$yList <- lapply(samp4$yList, function(x) x + rnorm(length(x)))
rcov4 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
brcov4 <- BinRawCov(rcov4)

h4 <- gcvlwls2d(sort(unique(unlist(samp4$tList))), kern='span', rcov=rcov4)
h4bin <- gcvlwls2d(sort(unique(unlist(samp4$tList))), kern='span', rcov=brcov4)
# fit4 <- lwls2d(h4$h, kern='epan', xin=rcov4$tpairn, yin=rcov4$cxxn, returnFit=TRUE)
# plot(fit4)

test_that('GCV on Binned rcov = rcov', {
  expect_equal(h3, h3bin)
  expect_equal(h4, h4bin)
})

# CV
set.seed(3)
system.time(h5 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), rcov=rcov3, kern='epan', CV='10fold'))
fit5 <- lwls2d(h5$h, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, returnFit=TRUE)
plot(fit5)

lcv(fit3)
tmp <- caret::createFolds(rcov3$cxxn, k=10)

# Consistency test: slow
set.seed(1)
pts <- seq(0, 1, by=0.01)
outPts <- seq(0, 1, by=0.1)
samp5 <- wiener(500, pts, sparsify=20)
rcov5 <- GetRawCov(samp5$yList, samp5$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
brcov5 <- BinRawCov(rcov5)
system.time(h5bin <- gcvlwls2d(sort(unique(unlist(samp5$tList))), kern='epan', rcov=brcov5))
# system.time(h5 <- gcvlwls2d(sort(unique(unlist(samp5$tList))), kern='epan', rcov=rcov5))
tmp <- lwls2d(h5bin$h, 'epan', brcov5$tPairs, brcov5$meanVals, brcov5$count, xout1=outPts, xout2=outPts)
test_that('gcvlwls2d smooth it well', {
  expect_equal(diag(tmp), outPts, tolerance=0.1)
})

# Consistency test: Dense
set.seed(1)
pts <- seq(0, 1, by=0.01)
outPts <- seq(0, 1, by=0.1)
samp6 <- wiener(10000, pts, sparsify=length(pts))
rcov6 <- GetRawCov(samp6$yList, samp6$tList, pts, rep(0, length(pts)), 'Dense', error=FALSE)
brcov6 <- BinRawCov(rcov6)
system.time(h6bin <- gcvlwls2d(sort(unique(unlist(samp6$tList))), kern='epan', rcov=brcov6))
tmp <- lwls2d(h6bin$h, 'epan', brcov6$tPairs, brcov6$meanVals, brcov6$count, xout1=outPts, xout2=outPts)
test_that('gcvlwls2d smooth it well', {
  expect_equal(diag(tmp), outPts, tolerance=0.02)
})
