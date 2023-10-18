
test_that('test the confidence level R', {
  set.seed(1000)
  # Test 1: scalar response ~ Functional + scalar predictors
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  lambdaX2 <- c(0.9,0.5)
  
  # training set
  n <- 50
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  Xi2 <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  sparseLt <- list(); sparseLy <- list()
  sparseLt2 <- list(); sparseLy2 <- list()
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLy <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[ind]
    
    denseLy <- lambdaX2[1]*Xi2[i,1]*phi1(t0,1) + lambdaX2[2]*Xi2[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt2[[i]] <- t0[ind]
    sparseLy2[[i]] <- denseLy[ind]
  }
  
  sparseX <- list(X1 = list(Ly=sparseLy,Lt=sparseLt), 
                  X2 = list(Ly=sparseLy2,Lt=sparseLt2), 
                  Z = rnorm(n,10,2))
  
  beta <- c(1, -1)
  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
  
  testthat::expect_error(FLMCI(Y=Y, X=sparseX, level=2, optnsListX=list(FVEthreshold=0.95), R = 100))
  testthat::expect_error(FLMCI(Y=Y, X=sparseX, level=-1, optnsListX=list(FVEthreshold=0.95), R = 100))
})


test_that('test the Bootstrap repetition R', {
  set.seed(1000)
  # Test 1: scalar response ~ Functional + scalar predictors
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  lambdaX2 <- c(0.9,0.5)
  
  # training set
  n <- 50
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  Xi2 <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  sparseLt <- list(); sparseLy <- list()
  sparseLt2 <- list(); sparseLy2 <- list()
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLy <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[ind]
    
    denseLy <- lambdaX2[1]*Xi2[i,1]*phi1(t0,1) + lambdaX2[2]*Xi2[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt2[[i]] <- t0[ind]
    sparseLy2[[i]] <- denseLy[ind]
  }
  
  sparseX <- list(X1 = list(Ly=sparseLy,Lt=sparseLt), 
                  X2 = list(Ly=sparseLy2,Lt=sparseLt2), 
                  Z = rnorm(n,10,2))
  
  beta <- c(1, -1)
  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
  
  testthat::expect_error(FLMCI(Y=Y, X=sparseX, level=2, optnsListX=list(FVEthreshold=0.95), R = 100.5))
})


test_that('scalar response ~ Functional + scalar predictors, CI output', {
  set.seed(1000)
  # Test 1: scalar response ~ Functional + scalar predictors
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  lambdaX2 <- c(0.9,0.5)
  
  # training set
  n <- 50
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  Xi2 <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  sparseLt <- list(); sparseLy <- list()
  sparseLt2 <- list(); sparseLy2 <- list()
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLy <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[ind]
    
    denseLy <- lambdaX2[1]*Xi2[i,1]*phi1(t0,1) + lambdaX2[2]*Xi2[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt2[[i]] <- t0[ind]
    sparseLy2[[i]] <- denseLy[ind]
  }
  
  sparseX <- list(X1 = list(Ly=sparseLy,Lt=sparseLt), 
                  X2 = list(Ly=sparseLy2,Lt=sparseLt2), 
                  Z = rnorm(n,10,2))
  
  beta <- c(1, -1)
  Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
  
  CI1 <- FLMCI(Y=Y,X=sparseX,optnsListX=list(FVEthreshold=0.95), R = 200)
  
  testthat::expect_lte(sum(CI1$CI_alpha$CI_upper<0), 1)
  testthat::expect_lte(sum(CI1$CI_alpha$CI_lower>0), 1)
  testthat::expect_lte(sum(CI1$CI_beta[[1]]$CI_lower>0), 5)
  testthat::expect_lte(sum(CI1$CI_beta[[1]]$CI_upper<0), 5)
  testthat::expect_lte(sum(CI1$CI_beta[[2]]$CI_lower>0), 5)
  testthat::expect_lte(sum(CI1$CI_beta[[2]]$CI_upper<0), 5)

})


test_that('functional response ~ Functional + scalar predictors', {
  set.seed(1000)
  ### functional covariate
  phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
  phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
  
  lambdaX <- c(1,0.7)
  lambdaX2 <- c(0.9,0.5)
  
  # training set
  n <- 50
  Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
  Xi2 <- matrix(rnorm(2*n),nrow=n,ncol=2)
  
  kxi = matrix(rnorm(2*n),nrow=n,ncol=2)
  
  sparseLt <- list(); sparseLy <- list()
  sparseLt2 <- list(); sparseLy2 <- list()
  YLt <- list(); YLy <- list()
  
  t0 <- seq(0,1,length.out=51)
  for (i in 1:n) {
    denseLy <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt[[i]] <- t0[ind]
    sparseLy[[i]] <- denseLy[ind]
    
    denseLy <- lambdaX2[1]*Xi2[i,1]*phi1(t0,1) + lambdaX2[2]*Xi2[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),3))
    sparseLt2[[i]] <- t0[ind]
    sparseLy2[[i]] <- denseLy[ind]
    
    denseLy <- lambdaX2[1]*kxi[i,1]*phi1(t0,1) + lambdaX2[2]*kxi[i,2]*phi1(t0,2)
    ind <- sort(sample(1:length(t0),5))
    YLt[[i]] <- t0[ind]
    YLy[[i]] <- denseLy[ind]
  }
  
  sparseX <- list(X1 = list(Ly=sparseLy,Lt=sparseLt), 
                  X2 = list(Ly=sparseLy2,Lt=sparseLt2), 
                  Z = rnorm(n,10,2))
  Y = list(Ly = YLy, Lt = YLt)
  flm_est = FLM1(Y=Y,X=sparseX,optnsListX=list(FVEthreshold=0.95))
  
  CI2 = FLMCI(Y=Y,X=sparseX,optnsListX=list(FVEthreshold=0.95), R = 200)
  
  testthat::expect_lte(sum(CI2$CI_alpha$CI_upper<0), 5)
  testthat::expect_lte(sum(CI2$CI_alpha$CI_lower>0), 5)
  testthat::expect_lte(sum(CI2$CI_beta[[1]]$CI_lower>0), 5)
  testthat::expect_lte(sum(CI2$CI_beta[[1]]$CI_upper<0), 5)
  testthat::expect_lte(sum(CI2$CI_beta[[2]]$CI_lower>0), 5)
  testthat::expect_lte(sum(CI2$CI_beta[[2]]$CI_upper<0), 5)
lao})


