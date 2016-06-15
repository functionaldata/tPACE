library(testthat)
# devtools::load_all()

set.seed(1)
n <- 300
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- sampWiener + matrix(rnorm(n, sd=10), n, length(pts))
sampWiener <- Sparsify(sampWiener, pts, 1:5)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', kernel='epan', 
                 methodBwCov='GCV'))

test_that('CreatePathPlot works', {
  CreatePathPlot(res)
  CreatePathPlot(res, subset=seq_len(n) %% 5 == 0, k=5, inputData=list(Lt=sampWiener$Lt, Ly=sampWiener$Ly), main='123', xlab='T')
})

test_that('User defined colors work', {
  showInd <- 11:13
  shown <- length(showInd)
  CreatePathPlot(res, showInd)
  CreatePathPlot(res, showInd, col=c('blue', 'cyan', 'grey'))
})
