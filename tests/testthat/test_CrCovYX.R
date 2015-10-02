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
