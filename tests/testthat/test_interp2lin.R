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
  BB = interp2lin(x, y, Z, xp, yp)
  expect_equal(sum(AA), sum(BB), tolerance = 2e-15) 

})



# These check out OK.
U = test_that("basic pracma::interp2 (extended) example gives same output ", {

  x <- linspace(-1, 4, 11)
  y <- linspace(-1, 4, 11)
  mgrid <- meshgrid(x, y)
  Z <- mgrid$X^2 + 2*mgrid$Y + 3
  xp <- yp <- linspace(-1, 4, 51)
  method <- "linear"
  AA = interp2(x, y, Z, xp, yp, method)
  BB = interp2lin( x, y, Z, xp, yp)
  expect_equal(sum(AA), sum(BB), tolerance = 2e-15)

})


 
