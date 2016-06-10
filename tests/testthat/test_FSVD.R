devtools::load_all()
library(testthat)
library(mvtnorm)

### Simulaiton setting in FSVD paper
# Generate X, Y functional samples with certain covariance structure
K = 3 # number of singular pairs
t = 10 # domain end point
obsGrid = seq(0, t, length.out = 51)
# singular functions
# X
phi1 = sqrt(2/t) * sin(2*pi*obsGrid/t)
phi2 = -sqrt(2/t) * cos(4*pi*obsGrid/t)
phi3 = -sqrt(2/t) * cos(2*pi*obsGrid/t)
Phi = cbind(phi1, phi2, phi3)
# Y
psi1 = phi3; psi2 = phi1; psi3 = phi2
Psi = cbind(psi1, psi2, psi3)
# singular components
CovscX = matrix(c(8,3,-2,3,4,1,-2,1,3), nrow=3)
CovscY = matrix(c(6,-2,1,-2,4.5,1.5,1,1.5,3.25), nrow=3)
CovscXY = diag(c(3,1.5,0.5))
CXY = CovscXY[1,1]*phi1%*%t(psi1) + CovscXY[2,2]*phi2%*%t(psi2) + 
  CovscXY[3,3]*phi3%*%t(psi3)
rho1 = CovscXY[1,1]/sqrt(CovscX[1,1]*CovscY[1,1])
rho2 = CovscXY[2,2]/sqrt(CovscX[2,2]*CovscY[2,2])
rho3 = CovscXY[3,3]/sqrt(CovscX[3,3]*CovscY[3,3])
# observed error sigma2
sigma2 = 0.05

### test for consistency in dense case
test_that('consistent estimates for dense case', {
  n = 10000 # sample size
  set.seed(1)
  singscore = rmvnorm(n=n, mean = rep(0,6), sigma = cbind(rbind(CovscX,CovscXY),rbind(CovscXY, CovscY)))
  zeta1 = singscore[,1] # X
  zeta2 = singscore[,2] # X
  zeta3 = singscore[,3] # X
  xi1 = singscore[,4] #Y
  xi2 = singscore[,5] #Y
  xi3 = singscore[,6] #Y

  Xmat = t(Phi %*% t(singscore[,1:3])) + rnorm(n = n*length(obsGrid),mean = 0,sd = sqrt(sigma2))
  Ymat = t(Psi %*% t(singscore[,4:6])) + rnorm(n = n*length(obsGrid),mean = 0,sd = sqrt(sigma2))
  Xins = MakeFPCAInputs(tVec = obsGrid, yVec = Xmat)
  Yins = MakeFPCAInputs(tVec = obsGrid, yVec = Ymat)
  
  Ly1 = Xins$Ly
  Lt1 = Xins$Lt
  Ly2 = Yins$Ly
  Lt2 = Yins$Lt
  fsvdobj = FSVD(Ly1, Lt1, Ly2, Lt2, SVDoptns = list(methodSelectK=3))
  expect_equal(fsvdobj$sValues, diag(CovscXY), tolerance = 0.01*diag(CovscXY))
  expect_equal(fsvdobj$canCorr, c(rho1, rho2, rho3), tolerance = 0.18*c(rho1, rho2, rho3))
  # test for high correlation between true and observed singular scores
  expect_more_than(abs(cor(fsvdobj$sScores1[,1], singscore[,1])), 0.99)
  expect_more_than(abs(cor(fsvdobj$sScores1[,2], singscore[,2])), 0.99)
  expect_more_than(abs(cor(fsvdobj$sScores1[,3], singscore[,3])), 0.99)
  expect_more_than(abs(cor(fsvdobj$sScores2[,1], singscore[,4])), 0.99)
  expect_more_than(abs(cor(fsvdobj$sScores2[,2], singscore[,5])), 0.99)
  expect_more_than(abs(cor(fsvdobj$sScores2[,3], singscore[,6])), 0.99)  
})

### test for consistency in dense case
test_that('consistent estimates for sparse case', {
  n = 2000 # sample size
  set.seed(123)
  singscore = rmvnorm(n=n, mean = rep(0,6), sigma = cbind(rbind(CovscX,CovscXY),rbind(CovscXY, CovscY)))
  zeta1 = singscore[,1] # X
  zeta2 = singscore[,2] # X
  zeta3 = singscore[,3] # X
  xi1 = singscore[,4] #Y
  xi2 = singscore[,5] #Y
  xi3 = singscore[,6] #Y
  
  Xmat = t(Phi %*% t(singscore[,1:3])) + rnorm(n = n*length(obsGrid),mean = 0,sd = sqrt(sigma2))
  Ymat = t(Psi %*% t(singscore[,4:6])) + rnorm(n = n*length(obsGrid),mean = 0,sd = sqrt(sigma2))
  XinSparse = Sparsify(samp = Xmat, pts = obsGrid, sparsity = sample(6:10,size=n,replace=TRUE))
  YinSparse = Sparsify(samp = Ymat, pts = obsGrid, sparsity = sample(6:10,size=n,replace=TRUE))
  
  Ly1 = XinSparse$Ly
  Lt1 = XinSparse$Lt
  Ly2 = YinSparse$Ly
  Lt2 = YinSparse$Lt
  fsvdobj1 = FSVD(Ly1, Lt1, Ly2, Lt2, SVDoptns = list(methodSelectK=3,bw1=1.2,bw2=1.2))
  expect_equal(fsvdobj1$sValues, diag(CovscXY), tolerance = 0.13*diag(CovscXY))
  expect_equal(fsvdobj1$canCorr[1], rho1, tolerance = 0.15*rho1)
  # test for high correlation between true and observed singular scores
  expect_more_than(abs(cor(fsvdobj1$sScores1[,1], singscore[,1])), 0.975)
  expect_more_than(abs(cor(fsvdobj1$sScores1[,2], singscore[,2])), 0.91)
  expect_more_than(abs(cor(fsvdobj1$sScores1[,3], singscore[,3])), 0.95)
  expect_more_than(abs(cor(fsvdobj1$sScores2[,1], singscore[,4])), 0.975)
  expect_more_than(abs(cor(fsvdobj1$sScores2[,2], singscore[,5])), 0.97)
  expect_more_than(abs(cor(fsvdobj1$sScores2[,3], singscore[,6])), 0.92)
})
