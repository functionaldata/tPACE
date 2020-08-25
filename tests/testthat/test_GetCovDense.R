# devtools::load_all()
library(testthat)
##options(error=recover)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

# test_that('GetCovDense fails when data is too sparse', {
#   
#   set.seed(1)
#   n <- 200
#   p <- 101
#   pts <- seq(0, 1, length.out=p)
#   sigma2 <- 0.1
#   mu <- pts
#   samp <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE) + 
#           rnorm(n * length(pts), sd=sqrt(sigma2))
#   samp[seq_len(n - 1), 1] <- NA # only subject 1 was observed at the first time point
#   expect_error(GetCovDense(samp, colMeans(samp), list(error=FALSE, dataType='DenseWithMV')), 
#                "Data is too sparse to be considered DenseWithMV. Remove sparse observations or specify dataType='Sparse' for FPCA")
# 
# })


test_that('GetCovDense with noise, get sigma2', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  samp <- Wiener(n, pts) + matrix(pts, n, p, byrow=TRUE) + 
          rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- MakeFPCAInputs(tVec=pts, yVec=samp)

  optnsNoerr <- SetOptions(tmp$Ly, tmp$Lt, list(error=FALSE, dataType='Dense'))
  optnsErr <- SetOptions(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense'))
  noerr <- GetCovDense(samp, colMeans(samp), optnsNoerr)
  err <- GetCovDense(samp, colMeans(samp), optnsErr)

  eigNoerr <- GetEigenAnalysisResults(noerr$smoothCov, pts, optnsNoerr)
  eigErr <- GetEigenAnalysisResults(err$smoothCov, pts, optnsErr)

  expect_equal(err$sigma2, sigma2, tolerance=1e-2)
  expect_equal(eigNoerr$fittedCov, eigNoerr$fittedCov)
})
