devtools::load_all(); library(testthat)

test_that('noisy dense data, default arguments: ie. cross-sectional mu/cov, use IN score', {
  
  expect_more_than( abs( cor( Ksi[,1], fvpaResults$fpcaObjY$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Ksi[,2], fvpaResults$fpcaObjY$xiEst[,2])), 0.975)
  expect_more_than( abs( cor( Zeta[,1], fvpaResults$fpcaObjR$xiEst[,1])), 0.975)
  expect_more_than( abs( cor( Zeta[,2], fvpaResults$fpcaObjR$xiEst[,2])), 0.94)
  
})
 
library(testthat)
test_that('FindNN is correct', {
  tPairs <- matrix(c(0, 0,
                     1, 1,
                     0, 1,
                     1, 0,
                     0, 0.2,
                     0.2, 0,
                     0.2, 0.2,
                     0.2, 1,
                     1, 0.2), byrow=TRUE, ncol=2)
  expect_equal(FindNN(tPairs), 0.8)
})



test_that('BWNN works for large sample', {
  library(fdapace)
  set.seed(1)
  n <- 100
  pts <- seq(0, 1, length.out=100)
  sparsity <- 2:5
  samp <- Wiener(n, pts, sparsity)
  bw <- BwNN(samp[[1]]) # Lt or Lt
  expect_true(bw['cov'] >= 0.01 && bw['cov'] < 0.1 && bw['mu'] <= bw['cov'])
})

