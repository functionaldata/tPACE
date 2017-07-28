# devtools::load_all()
# This is an additional line
library(testthat)

test_that('noisy dense data, default arguments', {
  
  N = 44;
  eps = 0.23;
  M = 41;
  set.seed(123) 
  Tfinal = 3
  me <- function(t) exp(-Tfinal*(((t/Tfinal^2)-0.5))^2);
  T = seq(0,Tfinal,length.out = M) 
  recondingTimesMat = matrix(nrow = N, ncol = M)
  yMat = matrix(nrow = N, ncol = M)
  
  for (i in 1:N){
    peak = runif(min = 0.25,max =  0.75,1)*Tfinal 
    recondingTimesMat[i,] = Tfinal* sort( unique(c( seq(0.0 , peak, length.out = round((M+1)/2)),
                                                    seq( peak, Tfinal, length.out = round((M+1)/2))) ))
    yMat[i,] = me(recondingTimesMat[i,])* rnorm(1, mean=4.0, sd=  eps)    + rnorm(M, mean=0.0, sd=  eps) 
  }
  
  Y = as.list(as.data.frame(t(yMat)))
  X = rep(list(T),N)
  
  sss =  WFDA(Ly = Y, Lt = X )
  
  expect_gt( abs(  cor( colMeans(sss$aligned), me(T*3) ) ), 0.996) 
  expect_gt( abs(  min(apply( sss$aligned, 1, function(u) cor( u, me(T*3) ) )) ), 0.925)
  
})

test_that('noiseless dense data - only phase variations, default arguments', {
  
  N = 44;
  eps = 0.000;
  M = 41;
  set.seed(123) 
  Tfinal = 3
  me <- function(t) exp(-Tfinal*(((t/Tfinal^2)-0.5))^2);
  T = seq(0,Tfinal,length.out = M) 
  recondingTimesMat = matrix(nrow = N, ncol = M)
  yMat = matrix(nrow = N, ncol = M)
  
  for (i in 1:N){
    peak = runif(min = 0.25,max =  0.75,1)*Tfinal 
    recondingTimesMat[i,] = Tfinal* sort( unique(c( seq(0.0 , peak, length.out = round((M+1)/2)),
                                                    seq( peak, Tfinal, length.out = round((M+1)/2))) ))
    yMat[i,] = me(recondingTimesMat[i,])* rnorm(1, mean=4.0, sd=  eps)    + rnorm(M, mean=0.0, sd=  eps) 
  }
  
  Y = as.list(as.data.frame(t(yMat)))
  X = rep(list(T),N)
  
  sss =  WFDA(Ly = Y, Lt = X )
  
  expect_gt( abs(  cor( colMeans(sss$aligned), me(T*3) ) ), 0.996) 
  expect_gt( abs(  min(apply( sss$aligned, 1, function(u) cor( u, me(T*3) ) )) ), 0.965)
  
})

test_that('R implemenation are at least as good as current MATLAB WFPCA implementation in computing the mean',{
  
  load(system.file('testdata', 'YmatFromWFPCAexample_rng123.RData', package='fdapace'))
 
  Y = as.list(as.data.frame(t(QQ)))
  N = nrow(QQ)
  T = seq(0,1,length.out = ncol(QQ))
  X = rep(list(T),N)
  sss =  WFDA(Ly = Y, Lt = X )
  
  # In MATLAB
  # >> sum((mean(aligned) - me(T)).^2)
  # ans =
  #   0.0213
  # >> sum(abs(mean(aligned) - me(T)).^1)
  # ans =
  #   0.4686
  
  meTRUE = exp(-10*(T-0.5)^2); # True mean
  expect_lt(sum(abs(meTRUE- colMeans(sss$aligned))), 0.4686) 
  expect_lt(sum((meTRUE- colMeans(sss$aligned))^2), 0.0213) 
  
})