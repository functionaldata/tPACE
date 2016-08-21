# devtools::load_all('../../RPACE/tPACE')

test_that('legendre basis works', {
  M <- 20000
  pts <- seq(0, 1, length.out=M)
  K <- 10

  tmp <- CreateBasis(K, pts, 'legendre01')
  expect_equal(crossprod(tmp) / M, diag(K), scale=1, tol=1e-3)

  K <- 1
  tmp <- CreateBasis(K, pts, 'legendre01')
  expect_equal(crossprod(tmp) / M, matrix(1), scale=1, tol=1e-3)
})
