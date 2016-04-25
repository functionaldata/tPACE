devtools::load_all()
library(testthat)

test_that('The cross-covariance of two constant processes is zero.',{
  
  Ly1= list( rep(2.1,7), rep(2.1,3),2.1 );
  Lt1 = list(1:7,1:3, 1);
  Ly2 = list( rep(1.1,7), rep(1.1,3),1.1); 
  Lt2 = list(1:7,1:3, 1);
  Ymu1 = rep(55,7);
  Ymu2 = rep(1.1,7);
  
  AA<-GetCrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2)
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
  
  AA<-GetCrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2, bw1=2, bw2=2)
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
  
  AA<-GetCrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2, bw1=2, bw2=2)
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
  
  AA <- GetCrCovYX(Ly1 = yTrueB, Ly2 =yTrueA)
  
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
  
  ySparseA = Sparsify(yTrueA, s, c(3:5))    
  ySparseB = Sparsify(yTrueB, s, c(3:5))     
  
  BB1 <- GetCrCovYX(Ly1 = ySparseA$Ly, Lt1 = ySparseA$Lt, Ly2 = ySparseB$Ly, Lt2 = ySparseB$Lt, 
                Ymu1 = rep(0,M), Ymu2 = rep(0,M), useGAM = TRUE  )
  
  BB2 <- GetCrCovYX(Ly1 = ySparseA$Ly, Lt1 = ySparseA$Lt, Ly2 = ySparseB$Ly, Lt2 = ySparseB$Lt, 
                 Ymu1 = rep(0,M), Ymu2 = rep(0,M), bw1=0.4, bw2=0.4  )
  
  sSmall = seq(0,10,length.out = 51)
  
  # we know that the covariance between ksi_1 and ksi_2 is three
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct1(sSmall))*3 - BB1$smoothedCC )),  0.02, tol=.01, scale=1 )
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct1(sSmall))*3 - BB2$smoothedCC )),  0.02, tol=.01, scale=1 )
   
})

test_that('The cross-covariance of two simple unrelated process is correct. Same readings lengths.',{
  
  set.seed(123)
  N = 1511;   
  M = 101;
  
  # Define the continuum
  s = seq(0,10,length.out = M)
  
  # Define the mean and 2 eigencomponents 
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5) 
  eigFunct2 <- function(s) +sin(2*s*pi/10) / sqrt(5) 
  
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale) 
  Ksi = t(t(chol(matrix(c(5,3,3,4),2))) %*% t(Ksi))
  
  # Create Y_true
  yTrueA = Ksi[,1] %*% t(matrix(eigFunct1(s), ncol=1)) 
  yTrueB = Ksi[,2] %*% t(matrix(eigFunct2(s), ncol=1)) 
  
  ySparseA = Sparsify(yTrueA, s, c(3:5))    
  ySparseB = Sparsify(yTrueB, s, c(3:5))     
  
  BB1 <- GetCrCovYX(Ly1 = ySparseA$Ly, Lt1 = ySparseA$Lt, Ly2 = ySparseB$Ly, Lt2 = ySparseB$Lt, 
                 Ymu1 = rep(0,M), Ymu2 = rep(0,M), useGAM = TRUE  )
  
  BB2 <- GetCrCovYX(Ly1 = ySparseA$Ly, Lt1 = ySparseA$Lt, Ly2 = ySparseB$Ly, Lt2 = ySparseB$Lt, 
                 Ymu1 = rep(0,M), Ymu2 = rep(0,M), bw1=0.4, bw2=0.4  )
  
  sSmall = seq(0,10,length.out = 51)
  
  # we know that the covariance between ksi_1 and ksi_2 is three
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct2(sSmall))*3 - BB1$smoothedCC )),  0.02, tol=.01, scale=1 )
  expect_equal( median(abs( eigFunct1(sSmall)%*%t(eigFunct2(sSmall))*3 - BB2$smoothedCC )),  0.02, tol=.01, scale=1 )

  # par(mfrow(1,3))
  # plot3D::persp3D(s,s, z= cov(yTrueA,yTrueB))
  # plot3D::persp3D(sSmall, sSmall, BB1$smoothedCC )
  # plot3D::persp3D(sSmall, sSmall, BB2$smoothedCC )
  
})

test_that('Dense Wiener process has cov(s,t) = min(s,t)', {
  set.seed(4)
  n <- 500
  nGridIn <- 51
  sparsity <- 1:5 # must have length > 1
  bw <- NA
  T <- matrix(seq(0, 1, length.out=nGridIn))

## Corr(X(t), Y(t)) = 1/2
  A <- Wiener(n, T)
  B <- Wiener(n, T)
  C <- Wiener(n, T)
  X <- A + B
  Y <- A + C

  tmp <- GetCrCovYX(bw, bw, Ly1=X, Ly2=Y)
  tmp1 <- GetCrCovYX(NULL, NULL, Ly1=X, Ly2=Y)
  expect_equal(diag(tmp$rawCC$rawCCov), as.numeric(T), tolerance=0.1)
  expect_equal(tmp, tmp1) # for dense data no smoothing is used.
})

test_that('Sparse Wiener process has cov(s,t) = min(s,t)', {
  set.seed(4)
  n <- 200
  nGridIn <- 51
  sparsity <- 1:5 # must have length > 1
  bw <- 0.2
  kern <- 'epan'
  T <- matrix(seq(0, 1, length.out=nGridIn))

## Corr(X(t), Y(t)) = 1/2
  A <- Wiener(n, T)
  B <- Wiener(n, T)
  C <- Wiener(n, T)# + matrix((1:nGridIn) , n, nGridIn, byrow=TRUE)
  X <- A + B
  Y <- A + C
  indEach <- lapply(1:n, function(x) sort(sample(nGridIn, sample(sparsity, 1))))
  tAll <- lapply(1:n, function(i) T[indEach[[i]]])
  Xsp <- lapply(1:n, function(i) X[i, indEach[[i]]])
  Ysp <- lapply(1:n, function(i) Y[i, indEach[[i]]])

  tmp <- GetCrCovYX(bw, bw, Xsp, tAll, rep(0, nGridIn), Ysp, tAll, rep(0, nGridIn))
  tmpGCV <- GetCrCovYX(NULL, NULL, Xsp, tAll, rep(0, nGridIn), Ysp, tAll, rep(0, nGridIn))
  expect_equal(diag(tmp$smoothedCC), as.numeric(T), tolerance=0.37)
  expect_equal(diag(tmpGCV$smoothedCC), as.numeric(T), tolerance=0.37)
})
