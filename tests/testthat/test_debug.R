devtools::load_all()
# library(tPACE)
# library(testthat)
# library(microbenchmark)
 
matXinYinWin <- function(xin, yin, win) {
  N <- nrow(xin)
  x1 <- factor(xin[, 1])
  x2 <- factor(xin[, 2])
  yMat <- matrix(0, nlevels(x1), nlevels(x2))
  wMat <- yMat

  for (i in seq_len(N)) {
    yMat[x1[i], x2[i]] <- yMat[x1[i], x2[i]] + yin[i]
    wMat[x1[i], x2[i]] <- wMat[x1[i], x2[i]] + 1
  }
  yMat <- yMat / wMat

  list(x1Grid = as.numeric(levels(x1)),
       x2Grid = as.numeric(levels(x2)),
       yMat = yMat, wMat = wMat)
}

## time points are on a finite grid
set.seed(1)
n <- 1e5
bw <- 0.2
xin <- round(matrix(runif(2 * n), n, 2), 2)
yin <- rnorm(n)
win <- rep(1, n)
xout1 <- xout2 <- seq(0, 1, length.out=5)
system.time(matList <- matXinYinWin(xin, yin, win))

tmp <- invisible(with(matList, RmatLwls2d(c(bw, bw), 'epan', yMat, wMat, x1Grid, x2Grid, xout1, xout2, FALSE, FALSE)))
