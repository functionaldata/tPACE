devtools::load_all()
options(error=recover)
library(testthat)

# GetRawCov
load('../data/dataForGetRawCov.RData')
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output
# rm(t, y)

set.seed(1)
pts <- seq(0, 1, by=0.01)
samp1 <- wiener(100, pts, sparsify=5:10)
rcov1 <- GetRawCov(samp1$yList, samp1$tList, pts, rep(0, length(pts)), 'Sparse', FALSE)
sm1 <- lwls2d(0.3, 'epan', rcov1$tpairn, rcov1$cxxn, as.numeric(rcov1$win), returnFit=TRUE)
plot(sm1)


# lwls2d
lwls2d(0.5, 'epan', rcov$tpairn, rcov$cxxn, as.numeric(rcov$win), xout=dat())

tmp <- lwls2d(2, 'epan', rcov$tpairn, rcov$cxxn, as.numeric(rcov$win), xout1=seq(3, 9, by=0.5), xout2=seq(3, 9, by=0.5), returnFit=TRUE)
tmp <- lwls2d(1, 'epan', rcov$tpairn, rcov$cxxn, as.numeric(rcov$win), returnFit=TRUE)
gcv(tmp)
plot(tmp)

# getGCVscores
test_that('getGCVscores interface works', {
    expect_equal(
        getGCVscores(2, kern='epan', xin=rcov$tpairn, yin=rcov$cxxn, win=as.numeric(rcov$win)), 
        gcv(lwls2d(2, 'epan', rcov$tpairn, rcov$cxxn, as.numeric(rcov$win), returnFit=TRUE))['gcv']
    )
})

# gcvlwls2d
gcvlwls2d(t, rcov=rcov, kern='epan', h0=1)
gcvlwls2d(t, rcov=rcov, kern='epan', h0=2)
gcvlwls2d(t, rcov=rcov, kern='epan')

h1 <- gcvlwls2d(samp1$tList, rcov=rcov1, kern='epan', h0=0.1)
fit1 <- lwls2d(h1$h, kern='epan', xin=rcov1$tpairn, yin=rcov1$cxxn, returnFit=TRUE)
plot(fit1)


# use wiener process
set.seed(1)
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(500, pts, sparsify=20)
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
system.time(h3 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=rcov3))
fit3 <- lwls2d(h3$h, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, returnFit=TRUE)
plot(fit3)

# Error=TRUE
samp4 <- samp3
samp4$yList <- lapply(samp4$yList, function(x) x + rnorm(length(x)))
rcov4 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
h4 <- gcvlwls2d(sort(unique(unlist(samp4$tList))), kern='span', rcov=rcov4)
fit4 <- lwls2d(h4$h, kern='epan', xin=rcov4$tpairn, yin=rcov4$cxxn, returnFit=TRUE)
plot(fit4)

# CV
set.seed(3)
system.time(h5 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), rcov=rcov3, kern='epan', CV='10fold'))
fit5 <- lwls2d(h5$h, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, returnFit=TRUE)
plot(fit5)

lcv(fit3)
tmp <- caret::createFolds(rcov3$cxxn, k=10)

# Binned rcov = rcov
brcov3 <- BinRawCov(rcov3)
system.time(h3 <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=rcov3))
system.time(h3bin <- gcvlwls2d(sort(unique(unlist(samp3$tList))), kern='epan', rcov=brcov3))
test_that('GCV on Binned rcov = rcov', {
  expect_equal(h3, h3bin)
})
