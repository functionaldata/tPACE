library(locfit)
data(ethanol)

fitRef <- locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='epan')
tmp <- lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE)

# the 2D smoother does not match matlab version because it only uses elliptical window (rather than rectangular as in matlab PACE)
test_that('The interface passes the arguments correctly', {
    expect_equal(fitted(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE, scale=TRUE), ethanol[, -1]), fitted(fitRef))
    expect_equal(fitted(lwls2d(0.5, kern='gauss', ethanol[, -1], ethanol[, 1], returnFit=TRUE, scale=TRUE), ethanol[, -1]), fitted(locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='gauss')))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], scale=TRUE), fitted(fitRef))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], xout=ethanol[1:5, 2:3], scale=TRUE), predict(fitRef, ethanol[1:5, 2:3]))
})

tmpDat <- ethanol[1:5, 2:3]
colnames(tmpDat) <- c('x1', 'x2')
tmp <- lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], xout=ethanol[1:5, 2:3], returnFit=TRUE, scale=TRUE)
tmp1 <- lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE, scale=TRUE)
gcv(tmp)
# gcv(tmp, ev=tmpDat)

