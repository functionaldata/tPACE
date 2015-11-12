devtools::load_all()
library(testthat)

test_that('The cross-covariance of two constant processes is zero.',{
  
  Ly1= list( rep(2.1,7), rep(2.1,3),2.1 );
  Lt1 = list(1:7,1:3, 1);
  Ly2 = list( rep(1.1,7), rep(1.1,3),1.1); 
  Lt2 = list(1:7,1:3, 1);
  Ymu1 = rep(55,7);
  Ymu2 = rep(1.1,7);
  
  AA<-CrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2)
  expect_equal( 0,  sum(AA$smoothedCC)  )
  
})

test_that('The cross-covariance of two unrelated process is close to zero.',{
  
  N  = 100
  set.seed(123)
  Ly1 = lapply(1:N, function(x) runif(7))
  Ly2 = lapply(1:N, function(x) runif(7)) 
  Lt1 = lapply(1:N, function(x) sort(c(0, runif(5),1)) )
  Lt2 = Lt1
  Ymu1 = rep(0.5, length(unique(unlist(Lt1))))
  Ymu2 = rep(0.1^9,  length(unique(unlist(Lt2))))
  
  AA<-CrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2, bw1=2, bw2=2)
  expect_equal( 0.0,  mean(AA$smoothedCC), tol=1e-3 )
})

test_that('The cross-covariance of two unrelated process is close to zero. Different readings lengths.',{
  
  N  = 100
  set.seed(123)
  Ly1 = lapply(1:N, function(x) runif(7))
  Ly2 = lapply(1:N, function(x) runif(4))
  Lt1 = lapply(1:N, function(x) sort(c(0, runif(5),1)) )
  Lt2 = lapply(1:N, function(x) sort(c(0, runif(2),1)) )
  Ymu1 = rep(0.5, length(unique(unlist(Lt1))))
  Ymu2 = rep(0.1^9,  length(unique(unlist(Lt2))))
  
  AA<-CrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2, bw1=2, bw2=2)
  expect_equal( 0.0,  mean(AA$smoothedCC), tol=1e-3 )
})

test_that('The cross-covariance of two simple related process is correct. Same readings lengths.',{
  
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
  yTrueA = Ksi[,1] %*% t(matrix(eigFunct1(s), ncol=1)) 
  yTrueB = Ksi[,2] %*% t(matrix(eigFunct1(s), ncol=1)) 
  
  AA <- CrCovYX(Ly1 = yTrueB, Ly2 =yTrueA)
  
  # we know that the covariance between ksi_1 and ksi_2 is three
  expect_equal( max(abs( eigFunct1(s)%*%t(eigFunct1(s))*3 - AA$rawCC$rawCCov )),  0.01, tol=.01, scale=1 )
})

test_that('The cross-covariance of two simple related process is correct. Same readings lengths.',{
  
  set.seed(123)
  N = 611;   
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
  yTrueA = Ksi[,1] %*% t(matrix(eigFunct1(s), ncol=1)) 
  yTrueB = Ksi[,2] %*% t(matrix(eigFunct1(s), ncol=1)) 
  
  ySparseA = sparsify(yTrueA, s, c(3:5))    
  ySparseB = sparsify(yTrueB, s, c(3:5))     
  
  BB1 <- CrCovYX(Ly1 = ySparseA$yList, Lt1 = ySparseA$tList, Ly2 = ySparseB$yList, Lt2 = ySparseB$tList, 
                Ymu1 = rep(0,M), Ymu2 = rep(0,M), fast = TRUE  )
  
  BB2 <- CrCovYX(Ly1 = ySparseA$yList, Lt1 = ySparseA$tList, Ly2 = ySparseB$yList, Lt2 = ySparseB$tList, 
                 Ymu1 = rep(0,M), Ymu2 = rep(0,M), bw1=0.4, bw2=0.4  )
  
  sSmall = seq(0,10,length.out = 51)
  
  # we know that the covariance between ksi_1 and ksi_2 is three
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct1(sSmall))*3 - BB1$smoothedCC )),  0.02, tol=.01, scale=1 )
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct1(sSmall))*3 - BB2$smoothedCC )),  0.02, tol=.01, scale=1 )
   
})