library(testthat)
library(locfit)
source('../lwls2d.R')
data(ethanol)

fitRef <- locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='epan')
tmp <- lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE)

testthat('The interface passes the arguments correctly', {
    expect_equal(fitted(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], returnFit=TRUE), ethanol[, -1]), fitted(fitRef))
    expect_equal(fitted(lwls2d(0.5, kern='gauss', ethanol[, -1], ethanol[, 1], returnFit=TRUE), ethanol[, -1]), fitted(locfit(NOx~lp(C, E, h=0.5, deg=1, scale=TRUE), data=ethanol, kern='gauss')))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1]), fitted(fitRef))
    expect_equal(lwls2d(0.5, kern='epan', ethanol[, -1], ethanol[, 1], xout=ethanol[1:5, 2:3]), predict(fitRef, ethanol[1:5, 2:3]))
})


