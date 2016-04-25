devtools::load_all()
library(testthat)

trueLam <- 4 / ((2 * (1:50) - 1 ) * pi) ^ 2

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.01)
samp3 <- Wiener(50, pts, sparsify=length(pts))
mu3 <- rep(0, length(pts))

# without error
p0 <- SetOptions(samp3$Ly, samp3$Lt, list(maxK=50, FVEthreshold=1, dataType='Sparse', error=FALSE, kernel='epan'))
noErrBin <- GetSmoothedCovarSurface(samp3$Ly, samp3$Lt, mu3, pts, regGrid, p0, useBinnedCov=TRUE)
tmp <- GetEigenAnalysisResults(noErrBin$smoothCov, regGrid, p0)

# consistency test
test_that('Eigenvalues are close', {
  expect_equal((abs(tmp$lam - trueLam[1:length(tmp$lam)]) / trueLam[1:length(tmp$lam)] )[1:3], trueLam[1:3], tolerance=0.2)
})

# TEst integrate to one.
innerProd <- apply(tmp$phi, 2, function(lam1) 
                   apply(tmp$phi, 2, function(lam2) 
                         pracma::trapz(noErrBin$outGrid, lam1 * lam2)))
test_that('Eigenfunctions are orthonormal', {
  expect_equal(diag(innerProd), rep(1, tmp$kChoosen))
  expect_equal(innerProd[row(innerProd) != col(innerProd)], rep(0, length(innerProd) - nrow(innerProd)), tolerance=0.010)
})


# # with error
# p1 <- SetOptions(samp3$Ly, samp3$Lt, CreateOptions(dataType='Sparse', error=TRUE, kernel='epan'))
# Err <- GetSmoothedCovarSurface(samp3$Ly, samp3$Lt, mu3, pts, regGrid, p1, useBinnedCov=FALSE)

