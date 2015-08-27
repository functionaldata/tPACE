library(testthat)
rmsdiff <- function(x1,x2){ return (sqrt(mean( (x1-x2)^2   )))}

test_that('FPCAreg for scalar case returns correct results for functional predictor',{

  # Set the random seed, the number of subjects (N)
  # and the number of measurements per subjects (M)
  set.seed(123)
  N = 200;   
  M = 100;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
 
  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)
 
  betaFunct <- function(s) s  - 7 *exp(-(s-7)^2)
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))
 
  # Create Y_true
  yTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

  # Create beta_Func
  betaFunc = c(4,2) %*%  t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2))
  
  # Create scalar dep. variable a
  a = c(); for (i in 1:N) { a[i] = rnorm(sd=0.1,1) + trapzRcpp(s, yTrue[i,] * betaFunc); }
   
  L3 = makePACEinputs(IDs = rep(1:N, each=M),tVec=rep(s,N), t(yTrue) )
  FPCAdense = FPCA(y = L3$Ly, t = L3$Lt)   
  FREGres <- FPCAregScalar(FPCAdense, depVar= a, bootStrap=TRUE)
  expect_equal( 0,  rmsdiff(betaFunc,FREGres$betaFunc) , tol= 1e-3  )
}
)

test_that('FPCAreg for scalar case returns correct results for functional predictor and scalar predictors',{
 
  set.seed(123)
  N = 200;
  M = 100;

  # Define the continuum
  s = seq(0,10,length.out = M)

  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) s  + 10*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) -sin(2*s*pi/10) / sqrt(5)

  betaFunct <- function(s) s  - 7 *exp(-(s-7)^2)

  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))

  # Create Y_true
  yTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

  # Create beta_Func
  betaFunc = c(4,3) %*%  t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2))

  # Create two additional scalar predictors
  x1 <- rt(N,1)
  x2 <- rt(N,2)

  # Create scalar dep. variable a
  a = c(); for (i in 1:N) { a[i] = rnorm(sd=0.001,1) + trapzRcpp(s, yTrue[i,] * betaFunc) + x1[i]*2 + x2[i]*3  }

  L3 = makePACEinputs(IDs = rep(1:N, each=M),tVec=rep(s,N), t(yTrue) )
  FPCAdense = FPCA(y = L3$Ly, t = L3$Lt)
  FREGres <- FPCAregScalar(FPCAdense, depVar= a, bootStrap=TRUE, extVar= data.frame(x1,x2))
 
  expect_equal( 0,  rmsdiff(betaFunc,FREGres$betaFunc) , tol= 1e-5  )

})


