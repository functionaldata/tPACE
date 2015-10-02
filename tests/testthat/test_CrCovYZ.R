devtools::load_all()
library(testthat)

test_that('The cross-covariance in the case of a dense matrix against a constant vector is zero and in the case of a steadily increasing matrix against a steadily increasing vector is stable',{
  set.seed(1)
  n <- 100
  dccObj <- CrCovYZ(bw=1, Z= rep(4,n), Ly=  matrix( runif(10*n), n))
  expect_equal( rep(0,10),  as.numeric(dccObj$rawCC$rawCCov) )
  dccObj <- CrCovYZ( Z= 1:n, Ly= matrix(1:(10*n),n))
  expect_equal( 0,  diff(range(dccObj$rawCC$rawCCov)) )
})

test_that('The cross-covariance in the case of sparse sample and constant vector is zero',{
  set.seed(1)
  # Make sparse sample
  yList <- list( runif(5),  c(1:3), c(2:4), c(4))
  tList <- list( c(1:5), c(1:3), c(1:3), 4)
  Z = rep(4,4)
  sccObj = CrCovYZ(bw=1, Z= Z, Ly=yList, Lt=tList, Ymu=rep(4,5))
  expect_equal( rep(0, sum(unlist(lapply(yList, length)))), sccObj$rawCC$rawCCov )
})

test_that('The cross-covariance in the case of a sparse sample that is steadily increasing and a vector that is steadily increasing is almost perfectly linear',{
  # Make sparse sample
  yList <- list(c(0:5),c(4.0000000000001), 5)
  tList <- list(c(0:5),c(4.0000000000001), 5)
  Z = c(2.5, 4.0000000000001, 5)
  sccObj = CrCovYZ( Z= Z, Ly=yList, Lt=tList, Ymu=rep(4.5,7))
  AA<- summary(lm( sccObj$smoothedCC ~  sort(unique(unlist(tList)))))
  expect_equal( AA$r.squared, 0.9998, tol=0.0001)
})


