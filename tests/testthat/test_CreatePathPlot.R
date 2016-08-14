library(testthat)
# devtools::load_all()

set.seed(1)
n <- 300
pts <- seq(0, 1, by=0.05)
sampWiener <- Wiener(n, pts)
sampWiener <- sampWiener + matrix(rnorm(n, sd=0.1), n, length(pts))
sampWiener <- Sparsify(sampWiener, pts, 1:5)
res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
            list(dataType='Sparse', kernel='epan'))
resDer <- FPCAder(res, list(method='DPC'))

test_that('CreatePathPlot works for FPCA object', {
  CreatePathPlot(res)
  CreatePathPlot(res, 1:10)
  CreatePathPlot(res, 1:20, showObs=FALSE)
  CreatePathPlot(res, 1:20, showMean=TRUE, showObs=FALSE)
  CreatePathPlot(res, 1:20, obsOnly=TRUE)
  CreatePathPlot(res, 1:20, obsOnly=TRUE, showObs=FALSE)
  CreatePathPlot(inputData=sampWiener, subset=1:20, obsOnly=TRUE)
  CreatePathPlot(res, subset=seq_len(n) %% 5 == 0, K=4, inputData=list(Lt=sampWiener$Lt, Ly=sampWiener$Ly), main='123', xlab='T')
})

test_that('CreatePathPlot works for FPCAder object', {
  CreatePathPlot(resDer)
  CreatePathPlot(resDer, 1:10)
  CreatePathPlot(resDer, 1:10, showMean=TRUE)
  CreatePathPlot(resDer, 1:20, showObs=TRUE)
  CreatePathPlot(resDer, 1:20, obsOnly=TRUE, showObs=FALSE)
})

test_that('User defined colors work', {
  showInd <- 11:13
  shown <- length(showInd)
  CreatePathPlot(res, showInd)
  CreatePathPlot(res, showInd, col=c('blue', 'cyan', 'grey'))
})
