devtools::load_all()
library(testthat)
##options(error=recover)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

test_that('GetCovDense with noise, get sigma2', {
  set.seed(1)
  n <- 200
  p <- 101
  pts <- seq(0, 1, length.out=p)
  sigma2 <- 0.1
  mu <- pts
  samp <- wiener(n, pts) + matrix(pts, n, p, byrow=TRUE) + 
          rnorm(n * length(pts), sd=sqrt(sigma2))
  tmp <- makeFPCAinputs(tVec=pts, yVec=samp)

  optnsNoerr <- SetOptions(tmp$Ly, tmp$Lt, list(error=FALSE, dataType='Dense'))
  optnsErr <- SetOptions(tmp$Ly, tmp$Lt, list(error=TRUE, dataType='Dense'))
  noerr <- GetCovDense(samp, colMeans(samp), optnsNoerr)
  err <- GetCovDense(samp, colMeans(samp), optnsErr)

  eigNoerr <- GetEigenAnalysisResults(noerr$smoothCov, pts, optnsNoerr)
  eigErr <- GetEigenAnalysisResults(err$smoothCov, pts, optnsErr)

  expect_equal(err$sigma2, sigma2, tolerance=1e-2)
  expect_equal(eigNoerr$fittedCov, eigNoerr$fittedCov)
})
