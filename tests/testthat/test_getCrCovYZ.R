# devtools::load_all()
library(testthat)

test_that('The cross-covariance in the case of a dense matrix against a constant vector is zero and in the case of a steadily increasing matrix against a steadily increasing vector is stable',{
  set.seed(1)
  n <- 100
  dccObj <- GetCrCovYZ(bw=1, Z= rep(4,n), Ly=  matrix( runif(10*n), n))
  expect_equal( rep(0,10),  as.numeric(dccObj$rawCC$rawCCov) )
  dccObj <- GetCrCovYZ( Z= 1:n, Ly= matrix(1:(10*n),n))
  expect_equal( 0,  diff(range(dccObj$rawCC$rawCCov)) )
})

test_that('The cross-covariance in the case of sparse sample and constant vector is zero',{
  set.seed(1)
  # Make sparse sample
  Ly <- list( runif(5),  c(1:3), c(2:4), c(4))
  Lt <- list( c(1:5), c(1:3), c(1:3), 4)
  Z = rep(4,4)
  sccObj = GetCrCovYZ(bw=1, Z= Z, Ly=Ly, Lt=Lt, Ymu=rep(4,5))
  expect_equal( rep(0, sum(unlist(lapply(Ly, length)))), sccObj$rawCC$rawCCov )
})

test_that('The cross-covariance in the case of a sparse sample that is steadily increasing and a vector that is steadily increasing is almost perfectly linear',{
  # Make sparse sample
  Ly <- list(c(0:5),c(4.0000000000001), 5)
  Lt <- list(c(0:5),c(4.0000000000001), 5)
  Z = c(2.5, 4.0000000000001, 5)
  sccObj = GetCrCovYZ( Z= Z, Ly=Ly, Lt=Lt, Ymu=rep(4.5,7))
  AA<- summary(lm( sccObj$smoothedCC ~  sort(unique(unlist(Lt)))))
  expect_equal( AA$r.squared, 0.9998, tol=0.001)
})

test_that('The cross-covariance in the case of dense sample and a random variable with known variance',{
  
  set.seed(123)
  N = 1111;   
  M = 101;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
  
  # Define the mean and 2 eigencomponents 
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5) 
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale) 
  Ksi = t(t(chol(matrix(c(5,3,3,4),2))) %*% t(Ksi))
  
  # Create Y_true
  yTrue = Ksi[,1] %*% t(matrix(eigFunct1(s), ncol=1)) 
  sccObj = GetCrCovYZ(Z= Ksi[,2], Ly=yTrue )
  
  # we know that the covariance between ksi_1 and z is three
  expect_equal( max( abs( eigFunct1(s)*3 - sccObj$rawCC$rawCCov)),  0.03, tol=.01, scale=1 )
}) 

test_that('The cross-covariance in the case of sparse sample and a random variable with known variance and the GCV bandwidth choice.',{
  
  set.seed(123)
  N = 3000;   
  M = 101;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
  
  # Define the mean and 2 eigencomponents 
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5) 
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale) 
  Ksi = t(t(chol(matrix(c(5,3,3,4),2))) %*% t(Ksi))
  
  # Create Y_true
  yTrue = Ksi[,1] %*% t(matrix(eigFunct1(s), ncol=1)) 
  ySparse = Sparsify(yTrue, s, c(3:9))    
  
  # Use GCV to pick the bandwidth
  sccObj = GetCrCovYZ( Z= Ksi[,2],Ly=ySparse$Ly, Lt=ySparse$Lt, Ymu = rep(0,M)  )
  
  # Uncomment to visually check coherence.  
  # plot(s, sccObj$smoothedCC)
  # lines(s, 3* eigFunct1(s))
  
  # we know that the covariance between ksi_1 and z is three
  expect_equal(  mean(abs( 3*eigFunct1(s) - sccObj$smoothedCC)), 0.035, tol=.1, scale=1 )
  
  # check that the relevant GCV scores are worse
  sccObjDOUBLE = GetCrCovYZ( bw = sccObj$bw*2,  Z= Ksi[,2],Ly=ySparse$Ly, Lt=ySparse$Lt, Ymu = rep(0,M)  )
  sccObjHALF = GetCrCovYZ( bw = sccObj$bw*0.5,  Z= Ksi[,2],Ly=ySparse$Ly, Lt=ySparse$Lt, Ymu = rep(0,M)  )

  expect_equal( min(c( sccObj$score, sccObjDOUBLE$score, sccObjHALF$score) ) , sccObj$score )
})

test_that('Dense Wiener process has cov(int X(s) ds, X(t)) = int min(s,t) ds', {
  set.seed(4)
  n <- 2000
  nGridIn <- 51
  T <- matrix(seq(0, 1, length.out=nGridIn))

## Corr(\int X(s) ds, X(t)) = \int min(s,t) ds
  covTrue <- rowMeans(outer(as.numeric(T), as.numeric(T), pmin))

  A <- Wiener(n, T)
  B <- Wiener(n, T)
  X <- A + B
  Z <- rowMeans(A)

  tmp <- GetCrCovYZ(NULL, Z, NULL, X, NULL, NULL, T)
  expect_equal(as.numeric(tmp$rawCC$rawCCov), covTrue, tolerance=0.1)
})

test_that('Sparse Wiener process has cov(int X(s) ds, X(t)) = int min(s,t) ds', {
  set.seed(4)
  n <- 2000
  nGridIn <- 51
  sparsity <- 1:5 # must have length > 1
  bw <- 0.2
  kern <- 'epan'
  T <- matrix(seq(0, 1, length.out=nGridIn))

## Corr(\int X(s) ds, X(t)) = \int min(s,t) ds
  covTrue <- rowMeans(outer(as.numeric(T), as.numeric(T), pmin))

  A <- Wiener(n, T)
  B <- Wiener(n, T)
  X <- A + B
  indEach <- lapply(1:n, function(x) sort(sample(nGridIn, sample(sparsity, 1))))
  tAll <- lapply(1:n, function(i) T[indEach[[i]]])
  Xsp <- lapply(1:n, function(i) X[i, indEach[[i]]])
  Z <- rowMeans(A)

  tmp <- GetCrCovYZ(bw, Z, NULL, Xsp, tAll, rep(0, nGridIn), T)
  expect_equal(tmp$smoothedCC, covTrue, tolerance=0.1)
})
