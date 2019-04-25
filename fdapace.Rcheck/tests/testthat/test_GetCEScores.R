# devtools::load_all()
#options(error=recover)
library(testthat)

# test GetIndCEScores
obsGrid <- seq(0, 1, by=0.1)
yVec <- c(1, -1) 
muVec <- c(0, 0) 
lamVec <- c(6, 1) 
phiMat <- diag(1, 2)
Sigma_Yi <- diag(10, 2)
test_that('GetIndCEScores works', {
  expect_equal(GetIndCEScores(yVec, muVec, lamVec, phiMat, Sigma_Yi)[1:2], 
               list(xiEst=matrix(c(0.6, -0.1)), xiVar=diag(c(2.4, 0.9))))
  expect_equal(GetIndCEScores(c(), muVec, lamVec, phiMat, Sigma_Yi)[1:2], # integer(0) is not really an empty vector, 'c()' is more obvious
               list(xiEst=matrix(NA, length(lamVec)), xiVar=matrix(NA, length(lamVec), length(lamVec))))
  expect_equal(GetIndCEScores(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd=1)$fittedY, matrix(0))
})

GetIndCEScores(yVec[1], muVec[1], lamVec, phiMat[1, , drop=FALSE], Sigma_Yi[1, 1, drop=FALSE], newyInd=1)

# Set up
set.seed(1)
n <- 20
pts <- seq(0, 1, by=0.05)
truncPts <- seq(0.2, 0.8, 0.05)
regGrid <- seq(0, 1, by=0.01)
sigma2 <- 0.4
samp1 <- Wiener(n, pts) + rnorm(length(pts) * n, sd=sigma2)
samp1 <- Sparsify(samp1, pts, 2:7)
mu1 <- rep(0, length(pts))
pNoTrunc <- SetOptions(samp1$Ly, samp1$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', userBwCov = 0.2))
smc1 <- GetSmoothedCovarSurface(samp1$Ly, samp1$Lt, mu1, pts, regGrid, pNoTrunc)
eig1 <- GetEigenAnalysisResults(smc1$smoothCov, regGrid, pNoTrunc)

# test GetMuPhiSig
# no truncation
phiObs <- ConvertSupport(regGrid, pts, phi=eig1$phi)
CovObs <- ConvertSupport(regGrid, pts, Cov=eig1$fittedCov)
notruncSamp <- TruncateObs(samp1$Ly, samp1$Lt, pts)
expect_equal(sapply(notruncSamp$Ly, length), sapply(samp1$Ly, length))
tmp <- GetMuPhiSig(notruncSamp$Lt, pts, mu1, phiObs, CovObs + diag(smc1$sigma2, length(pts)))
expect_equal(sapply(tmp, function(x) length(x$muVec)), sapply(samp1$Ly, length))


# truncation
test_that('Observations with length 0 produces NA in the xiEst, xiVar, and fittedY', {
  set.seed(1)
  n <- 20
  pts <- signif(seq(0, 1, by=0.05), 4)
  truncPts <- signif(seq(0.2, 0.8, 0.05), 4)
  regGrid <- seq(0, 1, by=0.01)
  sigma2 <- 0.4
  samp1 <- Wiener(n, pts) + rnorm(length(pts) * n, sd=sigma2)
  samp1 <- Sparsify(samp1, pts, 2:7)
  mu1 <- rep(0, length(pts))
  pNoTrunc <- SetOptions(samp1$Ly, samp1$Lt, list(dataType='Sparse', error=TRUE, kernel='epan', userBwCov = 0.2))
  smc1 <- GetSmoothedCovarSurface(samp1$Ly, samp1$Lt, mu1, pts, regGrid, pNoTrunc)
  eig1 <- GetEigenAnalysisResults(smc1$smoothCov, regGrid, pNoTrunc)
  phiObs <- ConvertSupport(regGrid, truncPts, phi=eig1$phi)
  CovObs <- ConvertSupport(regGrid, truncPts, Cov=eig1$fittedCov)
  truncSamp <- TruncateObs(samp1$Ly, samp1$Lt, truncPts)
  tmp <- GetMuPhiSig(truncSamp$Lt, truncPts, mu1[1:length(truncPts)], phiObs, CovObs + diag(smc1$sigma2, length(truncPts)))
  expect_equal(sapply(tmp, function(x) length(x$muVec)), sapply(truncSamp$Ly, length))
  tmp1 <- GetCEScores(truncSamp$Ly[17:20], truncSamp$Lt[17:20], list(verbose=TRUE), rep(0, length(truncPts)), truncPts, CovObs, eig1$lambda, phiObs, smc1$sigma2)
  
  expect_equal(tmp1[[1, 1]], matrix(NA, length(eig1$lambda)))
  expect_equal(tmp1[[2, 1]], matrix(NA, length(eig1$lambda), length(eig1$lambda)))
  expect_equal(tmp1[[3, 1]], matrix(NA, 0, 0))
}) 

# Test GetCEScores: compare to Matlab
test_that('GetCEScores for sparse case matches Matlab', {
  y <- list(c(1, 2), 4, c(0, 2, 3))
  t <- list(c(1.5, 2.5), 2, c(1, 1.5, 2.5))
  # mu <- rep(0, 7)
  mu <- seq(0, 3, length.out=7)
  obsGrid <- seq(0, 3, length.out=7)
  pts <- seq(0, 1, length.out=7)
  phi <- cbind(sin(2 * pi * pts), cos(2 * pi * pts))
  lambda <- c(6, 1)
  sigma2 <- 0.4
  fittedCov <- phi %*% diag(lambda) %*% t(phi)
  tmp <- GetCEScores(y, t, list(), mu, obsGrid, fittedCov, lambda, phi, sigma2)
# xiEst are the same
  expect_equal(t(do.call(cbind, tmp[1, ])), 
               matrix(c(0.709244942754497,  0.337643678160920,
                        -2.017923270954032, -0.194174757281554,
                        -1.011227267892009, -0.329341317365270), 3, 2,
                      byrow=TRUE))
# xiVar are the same
  expect_equal(unlist(tmp[2, ]), c(0.568965517241379, 0.149314724790420,
                                   0.149314724790420, 0.281609195402299,
                                   0.757281553398058, -0.504480817738509,
                                   -0.504480817738509, 0.951456310679612,
                                   0.341317365269461, 0.155573425829540,
                                   0.155573425829540, 0.281437125748503))
# fitted Y_i are the same 
  expect_equal(unlist(tmp[3, ]), c(1.162356321839080,  2.054597701149425, 3.844660194174757, 0.288922155688622, 1.829341317365270, 3.211077844311377))
})

# # Matlab code:
# y = cell(1, 3); 
# t = y;
# y{1} = [1, 2]; y{2} = [4]; y{3} = [0, 2, 3]; 
# t{1} = [1.5, 2.5]; t{2} = [2]; t{3} = [1, 1.5, 2.5];
# # mu = zeros(1, 7); 
# mu = linspace(0, 3, 7); 
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
# [xi_est, xi_var, ypred] = getScores1(y, t, mu, phi, lambda, sigma, sigma_new, noeig, error, method, shrink, out1, regular, rho)



# no error case. This observation has a singular Sigma_Yi, which should not work in theory.




# truncation
test_that('Noiseless example', {
  set.seed(1)
  n <- 100
  M <- 50
  K <- 5
  samp <- MakeGPFunctionalData(n, M, K=K, lambda=2^(-seq(1, K)))
  # matplot(t(samp$Y[1:10, ])) 
  sampList <- MakeFPCAInputs(tVec=samp$pts, yVec=samp$Y)
  resIN <- FPCA(sampList$Ly, sampList$Lt, list(methodXi='IN', lean=TRUE, error=FALSE, rho='no'))
  resCE <- FPCA(sampList$Ly, sampList$Lt, list(methodXi='CE', lean=TRUE, error=FALSE, rho='no'))
  
  expect_gt( abs( cor( samp$xi[,1] , resIN$xiEst[,1])),0.97)
  expect_gt( abs(cor( samp$xi[,2] , resIN$xiEst[,2])),0.97)
  expect_gt( abs(cor( samp$xi[,3] , resIN$xiEst[,3])),0.97)
  expect_gt( abs(cor( samp$xi[,4] , resIN$xiEst[,4])), 0.97)
  expect_gt( abs(cor( samp$xi[,5] , resIN$xiEst[,5])), 0.94)
   
  expect_gt( abs( cor( samp$xi[,1] , resCE$xiEst[,1])),0.97)
  expect_gt( abs(cor( samp$xi[,2] , resCE$xiEst[,2])),0.97)
  expect_gt( abs(cor( samp$xi[,3] , resCE$xiEst[,3])),0.97)
  expect_gt( abs(cor( samp$xi[,4] , resCE$xiEst[,4])), 0.97)
  expect_gt( abs(cor( samp$xi[,5] , resCE$xiEst[,5])), 0.94)
 
}) 

