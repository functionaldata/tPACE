if (Sys.getenv('TRAVIS') != 'true') {# Do not run on travis since this is slow
# devtools::load_all()
library(testthat)

# test for consistency in dense case
test_that('consistent estimates for dense case using IN', {
  
  set.seed(1)
  n <- 30
  pts <- seq(0, 1, length.out = 33)
  sampWiener <- Wiener(n, pts)
  sampWiener <- Sparsify(sampWiener, pts, 33)
  resA <- FPCA(sampWiener$Ly, sampWiener$Lt)
  testK = 4;
  
  AAA <- predict(resA, newLy = sampWiener$Ly,newLt =  sampWiener$Lt, K = testK, xiMethod = 'IN') 
  
  expect_lt( max(  abs(colMeans( abs(AAA$scores) - abs(resA$xiEst[,1:testK]) ) ) ), .Machine$double.eps)  
  expect_gt( min(  abs(diag(cor( AAA$scores , resA$xiEst[,1:testK])) ) ), 0.99999)  
  
})

# test for consistency in sparse case
test_that('consistent estimates for dense case using CE', {
  
  set.seed(1)
  n <- 30
  pts <- seq(0, 1, length.out = 33)
  sampWiener <- Wiener(n, pts)
  sampWiener <- Sparsify(sampWiener, pts, 33)
  resB <- FPCA(sampWiener$Ly, sampWiener$Lt)
  testK = 4;
  BBB <- predict(resB, newLy = sampWiener$Ly,newLt =  sampWiener$Lt, K = testK, xiMethod = 'CE')

  expect_lt( max( abs(colMeans( abs(BBB$scores) - abs(resB$xiEst[,1:testK])) ) ), sd(abs(resB$xiEst[,1:testK])) / 50)
  expect_gt( min( abs(diag(cor(BBB$scores , resB$xiEst[,1:testK])) ) ), 0.99)  
  
})

# test for consistency in sparse case
test_that('consistent estimates for sparse case using CE', {
  
  set.seed(1)
  n <- 70
  pts <- seq(0, 1, by=0.0025)
  sampWiener <- Wiener(n, pts)
  sampWiener <- Sparsify(sampWiener, pts, 4)
  resC <- FPCA(sampWiener$Ly, sampWiener$Lt)
  testK = 4
  
  CCC <- predict(resC, newLy = sampWiener$Ly, newLt =  sampWiener$Lt, K = testK)
  
  expect_lt( max( abs(colMeans(CCC$scores - resC$xiEst[,1:testK]) ) - 1.96 * apply(CCC$scores - resC$xiEst[,1:testK],2,sd) ), .Machine$double.eps)  
  expect_gt( min( abs(diag(cor(CCC$scores , resC$xiEst[,1:testK])) ) ), 0.9999)  
  
})

# test for consistency when using dense data with a sparse object
test_that('consistent estimates for dense data with sparse FPCA object', {
  
  set.seed(1)
  n <- 1001 # This object has to be pretty informed...
  pts <- seq(0, 1, length.out = 33)
  sampWiener <- Wiener(n, pts)
  
  set.seed(2)
  sampWienerS <- Sparsify(sampWiener, pts, 9)
  resS <- FPCA(sampWienerS$Ly, sampWienerS$Lt)
  set.seed(2)
  sampWienerD <- Sparsify(sampWiener, pts, 33)
  resD <- FPCA(sampWienerD$Ly, sampWienerD$Lt)
  testK = 3
  
  AAA2 <- predict(resS, newLy = sampWienerD$Ly,newLt =  sampWienerD$Lt, K = testK)
  
  expect_lt( max( abs(colMeans( abs(AAA2$scores) - abs(resD$xiEst[,1:testK])) ) - 1.96 * apply( abs(AAA2$scores) - abs( resD$xiEst[,1:testK]),2,sd) ), .Machine$double.eps)  
  expect_gt( min( abs(diag(cor(AAA2$scores , resD$xiEst[,1:testK])) ) ), 0.9)
  
})

# test for consistency when using sparse data with a dense object
test_that('consistent estimates for sparse data with dense FPCA object', {
  
  set.seed(1)
  n <- 51
  pts <- seq(0, 1, length.out = 99)
  sampWiener <- Wiener(n, pts)
  
  set.seed(2)
  sampWienerS <- Sparsify(sampWiener, pts, 11)
  set.seed(2)
  sampWienerD <- Sparsify(sampWiener, pts, 99)
  resD <- FPCA(sampWienerD$Ly, sampWienerD$Lt)
  testK = 4
  
  AAA3 <- predict(resD, newLy = sampWienerS$Ly,newLt =  sampWienerS$Lt, K = testK, xiMethod = 'CE')
  
  expect_lt( max( abs(colMeans( abs(AAA3$scores) - abs(resD$xiEst[,1:testK])) ) - 1.96 * apply( abs(AAA3$scores) - abs( resD$xiEst[,1:testK]),2,sd) ), .Machine$double.eps)  
  expect_gt( min( abs(diag(cor(AAA3$scores, resD$xiEst[,1:testK])) ) ), 0.90) # 
  
})

# test for consistency when using sparse data with a dense object
test_that('error when using Tr. Num. Int with sparse data FPCA object', {
  
  set.seed(1)
  n <- 51
  pts <- seq(0, 1, length.out = 99)
  sampWiener <- Wiener(n, pts)
  
  set.seed(2)
  sampWienerS <- Sparsify(sampWiener, pts, 2)
  resS2 <- FPCA(sampWienerS$Ly, sampWienerS$Lt)

  expect_error(  predict(resS2, newLy = sampWienerS$Ly,newLt =  sampWienerS$Lt, K = testK, xiMethod = 'IN'), 
                 message =  'Trapezoid Numerical intergration (IN) is invalid for sparse data.') # 
  
})

}
