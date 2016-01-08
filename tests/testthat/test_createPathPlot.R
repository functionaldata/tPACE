library(testthat)
devtools::load_all()

set.seed(1)
n <- 30
pts <- seq(0, 1, by=0.05)
sampWiener <- wiener(n, pts)
sampWiener <- sparsify(sampWiener, pts, 5:10)
res <- FPCA(sampWiener$yList, sampWiener$tList, 
            list(dataType='Sparse', error=FALSE, kernel='epan',
                 verbose=TRUE))

test_that('createPathPlot works', {
  createPathPlot(res)
  createPathPlot(res, subset=seq_len(n) %% 5 == 0, k=5, inputData=list(t=sampWiener$tList, y=sampWiener$yList), main='123', xlab='T')
})
