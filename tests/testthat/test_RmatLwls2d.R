devtools::load_all()
library(testthat)
library(microbenchmark)
 
matXinYinWin <- function(xin, yin, win) {
  N <- nrow(xin)
  x1 <- factor(xin[, 1])
  x2 <- factor(xin[, 2])
  yMat <- matrix(0, nlevels(x1), nlevels(x2))
  wMat <- yMat

  for (i in seq_len(N)) {
    yMat[x1[i], x2[i]] <- yMat[x1[i], x2[i]] + yin[i] * win[i]
    wMat[x1[i], x2[i]] <- wMat[x1[i], x2[i]] + win[i]
  }
  yMat <- yMat / wMat

  list(x1Grid = as.numeric(levels(x1)),
       x2Grid = as.numeric(levels(x2)),
       yMat = yMat, wMat = wMat)
}

## time points are unique
set.seed(1)
n <- 1e3
bw <- 0.2
xin <- round(matrix(runif(2 * n), n, 2), 12)
yin <- rnorm(n)
win <- rep(1, n)
xout1 <- xout2 <- seq(0, 1, length.out=5)
# system.time(tmp <- Rmullwlsk(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE))
system.time(matList <- matXinYinWin(xin, yin, win))

test_that('matXinYinWin does the right long-to-mat conversion', {
  sapply(seq_len(n), function(i) {
    x1Ind <- which.min(abs(matList[['x1Grid']] - xin[i, 1]))
    x2Ind <- which.min(abs(matList[['x2Grid']] - xin[i, 2]))
    expect_equal(yin[i], matList[['yMat']][x1Ind, x2Ind])
  })
})
tmp <- invisible(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)))
tmp1 <- RmullwlskCC(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE)

test_that('Mat 2D smoother is the same as the long 2D smoother, unique time points', {
  expect_equal(tmp, t(tmp1))
})

## time points are on a finite grid
set.seed(1)
n <- 1e4
bw <- 0.1
xin <- round(matrix(runif(2 * n), n, 2), 3)
yin <- rnorm(n)
win <- rep(1, n)
xout1 <- xout2 <- seq(0, 1, length.out=5)
system.time(matList <- matXinYinWin(xin, yin, win))

tmp <- invisible(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)))
tmp1 <- RmullwlskCC(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE)
test_that('Mat 2D smoother is the same as the long 2D smoother', {
  expect_equal(tmp, t(tmp1))
})

# microbenchmark(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)), times=100L)
# microbenchmark(RmullwlskCC(c(bw, bw), 'epan', t(xin), yin, win, xout1, xout2, FALSE), times=100L)


## Functional setting
set.seed(1)
bw <- 0.1
pts <- seq(0, 1, length.out=100L)
mu <- rep(0, length(pts))
a <- wiener(1000, pts, sparsify=5:10)
xout1 <- xout2 <- seq(0, 1, length.out=20)
rcov <- BinRawCov(GetRawCov(a$yList, a$tList, pts, mu, 'Sparse', 0))
matList <- matXinYinWin(rcov$tPairs, rcov$meanVals, rcov$count)
tmp <- with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE))
tmp1 <- with(rcov, RmullwlskCC(c(bw, bw), 'epan', t(tPairs), meanVals, count, xout1, xout2, FALSE))
expect_equal(tmp, t(tmp1))

# # speed test
microbenchmark(tmp <- with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)), times=10L)
microbenchmark(tmp1 <- with(rcov, RmullwlskCC(c(bw, bw), 'epan', t(tPairs), meanVals, count, xout1, xout2, FALSE)), times=10L)
