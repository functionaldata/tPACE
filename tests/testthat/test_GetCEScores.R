devtools::load_all()
options(error=recover)
library(testthat)

# test GetIndCEScores
obsGrid <- seq(0, 1, by=0.1)
yVec <- c(1, -1) 
muVec <- c(0, 0) 
lamVec <- c(6, 1) 
phiMat <- diag(1, 2)
Sigma_Yi <- diag(10, 2)
test_that('GetIndCEScores works', {
  expect_equal(GetIndCEScores(yVec, muVec, lamVec, phiMat, Sigma_Yi), 
               list(xiEst=matrix(c(0.6, -0.1)), xiVar=diag(c(2.4, 0.9))))
  expect_equal(GetIndCEScores(integer(0), muVec, lamVec, phiMat, Sigma_Yi), list(xiEst=NULL, xiVar=NULL))
})

# Set up
set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
truncPts <- seq(0.2, 0.8, 0.05)
regGrid <- seq(0, 1, by=0.01)
sigma2 <- 0.4
samp1 <- wiener(n, pts) + rnorm(length(pts) * n, sd=sigma2)
samp1 <- sparsify(samp1, pts, 2:7)
mu1 <- rep(0, length(pts))
pNoTrunc <- SetOptions(samp1$yList, samp1$tList, CreateOptions(dataType='Sparse', error=TRUE, kernel='epan'))
smc1 <- GetSmoothedCovarSurface(samp1$yList, samp1$tList, mu1, pts, regGrid, pNoTrunc)
eig1 <- GetEigenAnalysisResults(smc1$smoothCov, regGrid, pNoTrunc)

# test GetMuPhiSig
# no truncation
phiObs <- ConvertSupport(regGrid, pts, phi=eig1$phi)
CovObs <- ConvertSupport(regGrid, pts, Cov=eig1$fittedCov)
notruncSamp <- TruncateObs(samp1$yList, samp1$tList, pts)
expect_equal(sapply(notruncSamp$y, length), sapply(samp1$yList, length))
tmp <- GetMuPhiSig(notruncSamp$t, pts, mu1, phiObs, CovObs + diag(smc1$sigma2, length(pts)))
expect_equal(sapply(tmp, function(x) length(x$muVec)), sapply(samp1$yList, length))


# truncation
phiObs <- ConvertSupport(regGrid, truncPts, phi=eig1$phi)
CovObs <- ConvertSupport(regGrid, truncPts, Cov=eig1$fittedCov)

truncSamp <- TruncateObs(samp1$yList, samp1$tList, truncPts)
tmp <- GetMuPhiSig(truncSamp$t, truncPts, mu1[1:length(truncPts)], phiObs, CovObs + diag(smc1$sigma2, length(truncPts)))
expect_equal(sapply(tmp, function(x) length(x$muVec)), sapply(truncSamp$y, length))

# Test GetCEScores: compare to Matlab
y <- list(c(1, 2), 4, c(0, 2, 3))
t <- list(c(1.5, 2.5), 2, c(1, 1.5, 2.5))
mu <- rep(0, 7)
obsGrid <- seq(0, 3, length.out=7)
pts <- seq(0, 1, length.out=7)
phi <- cbind(sin(2 * pi * pts), cos(2 * pi * pts))
lambda <- c(6, 1)
sigma2 <- 0.4
fittedCov <- phi %*% diag(lambda) %*% t(phi)
tmp <- GetCEScores(y, t, list(), mu, obsGrid, fittedCov, lambda, phi, sigma2)
test_that('GetCEScores for sparse case matches Matlab', {
# xiEst are the same
  expect_equal(t(do.call(cbind, tmp[1, ])), 
               matrix(c(-2.463692959041937, -4.035846541908064,
                        -2.411388100357868, -0.646551724137931,
                        -0.388349514563107, -1.362275449101796), 3, 2))
# xiVar are the same
  expect_equal(unlist(tmp[2, ]), c(0.568965517241379, 0.149314724790420,
                                   0.149314724790420, 0.281609195402299,
                                   0.757281553398058, -0.504480817738509,
                                   -0.504480817738509, 0.951456310679612,
                                   0.341317365269461, 0.155573425829540,
                                   0.155573425829540, 0.281437125748503))
})

# # Matlab code:
# y = cell(1, 3); 
# t = y;
# y{1} = [1, 2]; y{2} = [4]; y{3} = [0, 2, 3]; 
# t{1} = [1.5, 2.5]; t{2} = [2]; t{3} = [1, 1.5, 2.5];
# mu = zeros(1, 7); 
# out1 = linspace(0, 3, 7); 
# pts = linspace(0, 1, 7); 
# phi = [sin(2 * pi * pts)', cos(2 * pi * pts)'];
# lambda = [6, 1]; 
# sigma = 0;
# sigma_new = 0.4;
# noeig = 2;
# error = 1;
# method = 'CE';
# shrink = 0;
# regular = 0;
# rho = 0;
# [xi_est, xi_var, ~] = getScores1(y, t, mu, phi, lambda, sigma, sigma_new, noeig, error, method, shrink, out1, regular, rho)



# no error case. This observation has a singular Sigma_Yi, which should not work in theory.
tmp1 <- GetCEScores(y[3], t[3], list(), mu, obsGrid, fittedCov, lambda, phi, 0)
