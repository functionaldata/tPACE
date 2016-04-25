devtools::load_all(); library(testthat)


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

