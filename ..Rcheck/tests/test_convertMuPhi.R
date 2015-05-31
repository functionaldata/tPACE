library(testthat)

test_that('1D interpolation and extrapolation', {
    expect_equal(mapX1d(1:3, 1:3, 0:5), 0:5)
    expect_equal(mapX1d(1:3, 1:3, c(1, 1, 2, 3)), c(1, 1, 2, 3))
    expect_equal(mapX1d(1:3, matrix(c(1:3, 2:4), ncol=2), 0:5), matrix(c(0:5, 1:6), ncol=2))
})



tList <- list(c(2, 4, 6), 2:4)
out1 <- 1:3
mu <- 1:3 
phi <- matrix(c(1:3, 5:7), ncol=2)

test_that('convertMuPhi', {
    expect_equal(convertMuPhi(tList, out1, mu, phi, 'Dense'), list(muSub=c(2, 4, 6), phiSub=matrix(c(2, 4, 6, 6, 8, 10), ncol=2)))
    expect_equal(convertMuPhi(tList, out1, mu, phi, 'Sparse'), 
        list(muSub=list(c(2, 4, 6), c(2, 3, 4)), phiSub=list(matrix(c(2, 4, 6, 6, 8, 10), ncol=2), matrix(c(2, 3, 4, 6, 7, 8), ncol=2)))
    )
})
