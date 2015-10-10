devtools::load_all()
library(testthat)
#options(error=recover)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

# set.seed(1)
# n <- 100
# pts <- seq(0, 1, by=0.05)
# samp3 <- wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
# samp3 <- sparsify(samp3, pts, 10)
# res <- FPCA(samp3$yList, samp3$tList, list(dataType='Sparse', useBins=TRUE))
# res$lambda / trueLam[1:length(res$lambda)]
# res$sigma2

# createCovPlot(res, 'Smoothed', FALSE)

test_that('Truncation works for FPCA Wiener process', {
  set.seed(1)
  n <- 100
  pts <- seq(0, 1, by=0.01)
  mu <- rep(0, length(pts))
  samp4 <- wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
  samp4 <- sparsify(samp4, pts, 10)
  samp5 <- samp4
  samp4$yList[[1]] <- samp4$tList[[1]] <- c(0, 1)
  pTrunc <- SetOptions(samp4$yList, samp4$tList, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0.03, 0.97), verbose=TRUE))
  pNoTrunc <- SetOptions(samp4$yList, samp4$tList, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))
  set.seed(1); res4 <- FPCA(samp4$yList, samp4$tList, pTrunc)
  set.seed(1); res5NoTrunc <- FPCA(samp5$yList, samp5$tList, pNoTrunc)
  set.seed(1); res4NoTrunc <- FPCA(samp4$yList, samp4$tList, pNoTrunc)

  expect_equal(res4[c('sigma2', 'bwMu', 'bwCov')], res4NoTrunc[c('sigma2', 'bwMu', 'bwCov')])
  expect_equal(min(max(abs(res4$xiEst[-1, 1] - res4NoTrunc$xiEst[-1, 1])), max(abs(res4$xiEst[-1, 1] + res4NoTrunc$xiEst[-1, 1]))), 0, tol=0.5)
  expect_equal(min(max(abs(res5NoTrunc$xiEst[-1, 1] - res4NoTrunc$xiEst[-1, 1])), max(abs(res5NoTrunc$xiEst[-1, 1] + res4NoTrunc$xiEst[-1, 1]))), 0, tol=0.05)
  expect_equal(nrow(res4$xiEst), nrow(res4NoTrunc$xiEst))
  expect_equal(length(res4$xiVar), length(res4NoTrunc$xiVar))
})

test_that('Missing values work for FPCA Wiener process', {
  set.seed(1)
  n <- 200
  pts <- seq(0, 1, by=0.01)
  mu <- rep(0, length(pts))
  samp4 <- wiener(n, pts) + rnorm(n * length(pts), sd=0.1)
  samp4 <- sparsify(samp4, pts, 10)
  samp4$yList[[1]] <- samp4$tList[[1]] <- c(0, 1)
  pNoTrunc <- SetOptions(samp4$yList, samp4$tList, list(dataType='Sparse', error=TRUE, kernel='epan', outPercent=c(0, 1), verbose=TRUE))

  samp4$yList[[2]] <- samp4$tList[[2]] <- c(0.1, 0.2, 0.5)
  set.seed(1); res4 <- FPCA(samp4$yList, samp4$tList, pNoTrunc)
  
  samp4$yList[[2]] <- c(NA, 0.2, 0.5)
  set.seed(1); res4NaN <- FPCA(samp4$yList, samp4$tList, pNoTrunc)

  samp4$yList[[2]] <-  c(0.1, 0.2, 0.5)
  samp4$tList[[2]] <-  c(NA, 0.2, 0.5)
  set.seed(1); res4NaN2 <- FPCA(samp4$yList, samp4$tList, pNoTrunc)

  expect_equal(res4[c('sigma2', 'bwMu', 'bwCov')], res4NaN[c('sigma2', 'bwMu', 'bwCov')])
  expect_equal(nrow(res4$xiEst), nrow(res4NaN$xiEst))
  expect_equal(length(res4$xiVar), length(res4NaN$xiVar))
  expect_equal(res4NaN$inputData$y[[2]], c(0.2,0.5))
  expect_equal(   res4$inputData$y[[2]], c(0.1,0.2,0.5))
  expect_equal(res4NaN2[c('sigma2', 'bwMu', 'lambda')], res4NaN[c('sigma2', 'bwMu', 'lambda')])

})
