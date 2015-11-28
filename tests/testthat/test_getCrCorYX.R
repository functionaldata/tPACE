library(testthat)

test_that('getCrCorYX works on a trivial example', {
  
  set.seed(123)
  A = matrix(rnorm(10*7), ncol=7); 
  B = matrix(rnorm(10*4), ncol=4);
  
  covA = cov(A)
  covB = cov(B)
  covAB = cov(A,B) 
  
  expect_equal(getCrCorYX(covAB, covA, covB),cor(A,B))  
  
  expect_equal(getCrCorYX(covAB, diag(covA), diag(covB)),cor(A,B))  
  
})

test_that('getCrCorYX works on a trivial example with a scalar', {
  
  set.seed(123)
  A = matrix(rnorm(101*7), ncol=7); 
  B = rnorm(101);
  
  covA = cov(A)
  covB = var(B)
  covAB = cov(A,B) 
  
  expect_equal(getCrCorYX(covAB, covA, covB),cor(A,B))  
  
  expect_equal(getCrCorYX(covAB, diag(covA), covB),cor(A,B))  
  
})

