options(error=recover)
library(locfit)
library(pracma)
source('../GetRawCov.R')
source('../mapX1d.R')
source('../IsRegular.R')

# GetRawCov
load('../data/dataForGetRawCov.RData')
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output
# rm(t, y)

source('../wiener.R')
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

debug(gcvlwls2d)
undebug(gcvlwls2d)

# use wiener process
set.seed(1)
pts <- seq(0, 1, by=0.05)
samp3 <- wiener(200, pts, sparsify=2:7)
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', FALSE)
h3 <- gcvlwls2d(samp3$tList, kern='epan', rcov=rcov3)
fit3 <- lwls2d(h3$h, kern='epan', xin=rcov3$tpairn, yin=rcov3$cxxn, returnFit=TRUE)
plot(fit3)