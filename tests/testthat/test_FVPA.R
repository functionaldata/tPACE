devtools::load_all(); library(testthat)

test_that('noisy dense data, default arguments: ie. cross-sectional mu/cov, use IN score', {
  
  # Define the continuum and basic constants
  set.seed(123); N = 333;   M = 101;
  s = seq(0,10,length.out = M)
  
  # Create S_true
  meanFunct <- function(s) s + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  Phi <- matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)
  ev <- c(5, 2)^2 
  Ksi = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(ev)) 
  
  Sij = Ksi %*% t(Phi) + t(matrix(rep(meanFunct(s),N), nrow=M))
  
  # Create Signs
  epsilonSign = ifelse( matrix(runif(N*M)-0.5, ncol = M) > 0, 1, -1)
  
  # Create residuals
  sigmaW = 0.05
  Wij = matrix(rnorm(N*M, sd = sigmaW/1), ncol = M)
  eigFunctA <- function(s) +cos(9*s*pi/10) / sqrt(5)
  eigFunctB <- function(s) -cos(1*s*pi)/sqrt(5) 
  Psi = matrix(c(eigFunctA(s),eigFunctB(s)), ncol=2)
  evR <- c(0.5, 0.3)^2 
  Zeta = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(evR))
  Vij = Zeta %*% t(Psi)
  
  # matplot(t(  exp(Vij + Wij)^0.5  ), t='l')
  Rij = epsilonSign * ( exp(Vij + Wij)^(0.5) )
  
  Xij = Sij +  Rij;
  
  tmp <- MakeFPCAInputs(tVec=s, yVec=Xij)
  
  fvpaResults = FVPA(y= tmp$Ly, t= tmp$Lt)
  
  expect_more_than( abs( cor( eigFunct1(s), fvpaResults$fpcaObjY$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunct2(s), fvpaResults$fpcaObjY$phi[,2])), 0.975)
  expect_more_than( abs( cor( eigFunctA(s), fvpaResults$fpcaObjR$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunctB(s), fvpaResults$fpcaObjR$phi[,2])), 0.975)
  
  expect_more_than( abs( cor( Ksi[,1], fvpaResults$fpcaObjY$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Ksi[,2], fvpaResults$fpcaObjY$xiEst[,2])), 0.975)
  expect_more_than( abs( cor( Zeta[,1], fvpaResults$fpcaObjR$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Zeta[,2], fvpaResults$fpcaObjR$xiEst[,2])), 0.94)
  
})

test_that('noisy dense data, CE scores and cross-sectional mu/cov', {
  
  # Define the continuum and basic constants
  set.seed(123); N = 333;   M = 101;
  s = seq(0,10,length.out = M)
  
  # Create S_true
  meanFunct <- function(s) s + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  Phi <- matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)
  ev <- c(5, 2)^2 
  Ksi = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(ev)) 
  
  Sij = Ksi %*% t(Phi) + t(matrix(rep(meanFunct(s),N), nrow=M))
  
  # Create Signs
  epsilonSign = ifelse( matrix(runif(N*M)-0.5, ncol = M) > 0, 1, -1)
  
  # Create residuals
  sigmaW = 0.05
  Wij = matrix(rnorm(N*M, sd = sigmaW/1), ncol = M)
  eigFunctA <- function(s) +cos(9*s*pi/10) / sqrt(5)
  eigFunctB <- function(s) -cos(1*s*pi)/sqrt(5) 
  Psi = matrix(c(eigFunctA(s),eigFunctB(s)), ncol=2)
  evR <- c(0.5, 0.3)^2 
  Zeta = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(evR))
  Vij = Zeta %*% t(Psi)
   
  Rij = epsilonSign * ( exp(Vij + Wij)^(0.5) )
  
  Xij = Sij +  Rij;
  
  tmp <- MakeFPCAInputs(tVec=s, yVec=Xij)
  
  fvpaResults = FVPA(y= tmp$Ly, t= tmp$Lt, optns = list(error=TRUE, FVEthreshold = 0.9, methodXi='CE' ))
  
  expect_more_than( abs( cor( eigFunct1(s), fvpaResults$fpcaObjY$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunct2(s), fvpaResults$fpcaObjY$phi[,2])), 0.975)
  expect_more_than( abs( cor( eigFunctA(s), fvpaResults$fpcaObjR$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunctB(s), fvpaResults$fpcaObjR$phi[,2])), 0.975)
  
  expect_more_than( abs( cor( Ksi[,1], fvpaResults$fpcaObjY$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Ksi[,2], fvpaResults$fpcaObjY$xiEst[,2])), 0.975)
  expect_more_than( abs( cor( Zeta[,1], fvpaResults$fpcaObjR$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Zeta[,2], fvpaResults$fpcaObjR$xiEst[,2])), 0.94)
  
})

test_that('noisy dense data, IN scores and smooth mu/cov', {
  
  # Define the continuum and basic constants
  set.seed(123); N = 333;   M = 101;
  s = seq(0,10,length.out = M)
  
  # Create S_true
  meanFunct <- function(s) s + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
  Phi <- matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)
  ev <- c(5, 2)^2 
  Ksi = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(ev)) 
  
  Sij = Ksi %*% t(Phi) + t(matrix(rep(meanFunct(s),N), nrow=M))
  
  # Create Signs
  epsilonSign = ifelse( matrix(runif(N*M)-0.5, ncol = M) > 0, 1, -1)
  
  # Create residuals
  sigmaW = 0.05
  Wij = matrix(rnorm(N*M, sd = sigmaW/1), ncol = M)
  eigFunctA <- function(s) +cos(9*s*pi/10) / sqrt(5)
  eigFunctB <- function(s) -cos(1*s*pi)/sqrt(5) 
  Psi = matrix(c(eigFunctA(s),eigFunctB(s)), ncol=2)
  evR <- c(0.5, 0.3)^2 
  Zeta = apply(matrix(rnorm(N*2), ncol=2), 2, scale) %*% diag(sqrt(evR))
  Vij = Zeta %*% t(Psi)
  
  Rij = epsilonSign * ( exp(Vij + Wij)^(0.5) )
  
  Xij = Sij +  Rij;
  
  tmp <- MakeFPCAInputs(tVec=s, yVec=Xij)
  
  fvpaResults = FVPA(y= tmp$Ly, t= tmp$Lt, optns = list(error=TRUE, FVEthreshold = 0.9, methodMuCovEst='smooth', userBwCov = 0.600 ))
  
  expect_more_than( abs( cor( eigFunct1(s), fvpaResults$fpcaObjY$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunct2(s), fvpaResults$fpcaObjY$phi[,2])), 0.975)
  expect_more_than( abs( cor( eigFunctA(s), fvpaResults$fpcaObjR$phi[,1])), 0.975)
  expect_more_than( abs( cor( eigFunctB(s), fvpaResults$fpcaObjR$phi[,2])), 0.975)
  
  expect_more_than( abs( cor( Ksi[,1], fvpaResults$fpcaObjY$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Ksi[,2], fvpaResults$fpcaObjY$xiEst[,2])), 0.975)
  expect_more_than( abs( cor( Zeta[,1], fvpaResults$fpcaObjR$xiEst[,1])), 0.97)
  expect_more_than( abs( cor( Zeta[,2], fvpaResults$fpcaObjR$xiEst[,2])), 0.91)
  
})
