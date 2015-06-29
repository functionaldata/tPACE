# setwd('misc/') 

library(testthat)

# These check out OK.
U = test_that("basic pracma::interp2 example gives same output ", { 

  x <- linspace(-1, 1, 11)
  y <- linspace(-1, 1, 11)
  mgrid <- meshgrid(x, y)
  Z <- mgrid$X^2 + mgrid$Y^2
  xp <- yp <- linspace(-1, 1, 101)
  method <- "linear" 
  AA = interp2(x, y, Z, xp, yp, method)
  BB = interp2lin(  mgrid$X, mgrid$Y, Z[1:length(Z)],  xp, yp)
  expect_equal(sum(AA), sum(BB), tolerance = 2e-14) 

})
 
