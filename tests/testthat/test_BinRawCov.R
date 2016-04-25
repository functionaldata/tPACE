devtools::load_all()
library(testthat)

# GMeanAndGCV

set.seed(1)
pts <- c(0, 1, 3:100) / 100
regGrid <- seq(0, 1, by=0.1)
samp3 <- Wiener(200, pts, sparsify=2:7)
p0 <- SetOptions(samp3$Ly, samp3$Lt, optns=list(dataType='Sparse', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))
rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, mu3, p0$dataType, error=p0$error)

brcov3 <- BinRawCov(rcov3)

test_that('BinRawCov works', {
  tPairs <- matrix(c(1, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2), ncol=2, byrow=TRUE)
  rcov <- list(tPairs=tPairs, cxxn=1:nrow(tPairs), error=FALSE)
  brcov <- BinRawCov(rcov)
  expect_equal(brcov$tPairs, matrix(c(1, 1, 2, 1, 1, 2, 2, 2), ncol=2, byrow=TRUE))
  expect_equal(brcov$meanVals, c(1, 2.5, 6, 4))
  
  rcov <- list(tPairs=tPairs, cxxn=1:nrow(tPairs), error=TRUE)
  brcov <- BinRawCov(rcov)
  expect_equal(brcov$tPairs, matrix(c(1, 1, 2, 1, 1, 2, 2, 2), ncol=2, byrow=TRUE))
  expect_equal(brcov$meanVals, c(1, 2.5, 6, 4))  
})