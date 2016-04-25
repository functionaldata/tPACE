devtools::load_all()
##options(error=recover)
library(testthat)

test_that('FitEigenValues works for binned rcov, error=TRUE', {
  set.seed(2)
  pts <- seq(0, 1, by=0.001)
  samp3 <- Wiener(10, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
  brcov3 <- BinRawCov(rcov3)
  phi <- cbind(sin((2 * 1 - 1) * pi * pts / 2), sin((2 * 2 - 1) * pi * pts / 2)) * sqrt(2)
  expect_equal(FitEigenValues(rcov3, pts, phi, 2), FitEigenValues(brcov3, pts, phi, 2))
})

test_that('FitEigenValues works for binned rcov, error=FALSE', {
  set.seed(2)
  pts <- seq(0, 1, by=0.001)
  samp3 <- Wiener(10, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  brcov3 <- BinRawCov(rcov3)
  phi <- cbind(sin((2 * 1 - 1) * pi * pts / 2), sin((2 * 2 - 1) * pi * pts / 2)) * sqrt(2)
  expect_equal(FitEigenValues(rcov3, pts, phi, 2), FitEigenValues(brcov3, pts, phi, 2))  
})

trueLambda <- 4 / (2 * (1:20) - 1)^2 / pi^2
test_that('FitEigenValues is consistent', {
  set.seed(2)
  pts <- seq(0, 1, by=0.05)
  samp3 <- Wiener(300, pts, sparsify=length(pts))
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=TRUE)
  phi <- cbind(sin((2 * 1 - 1) * pi * pts / 2), sin((2 * 2 - 1) * pi * pts / 2)) * sqrt(2)
  estLam <- FitEigenValues(rcov3, pts, phi, 2)
  expect_equal(estLam, trueLambda[1:length(estLam)], tolerance = 0.15)
})

# Test truncation
