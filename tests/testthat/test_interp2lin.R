# setwd('misc/') 
# devtools::load_all()

library(testthat)
library(pracma)
# These check out OK.
U1 = test_that("basic pracma::interp2 example gives same output ", { 

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
U2 = test_that("basic pracma::interp2 (extended) example gives same output ", {

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

# These check out OK.
U3 = test_that("basic pracma::interp2 (extended) example gives same output and out of sample", {

  x <- linspace(-1, 4, 11)
  y <- linspace(-1, 4, 11)
  mgrid <- meshgrid(x, y)
  Z <- mgrid$X^2 + 2*mgrid$Y + 3
  xp <- yp <- linspace(-1, 5, 51)
  method <- "linear"
  AA = interp2(x, y, Z, xp, yp, method)
  BB = interp2lin( x, y, Z, xp, yp)
  expect_equal(sum(AA, na.rm=TRUE), sum(BB, na.rm=TRUE), tolerance = 2e-15)
  expect_equal(is.na(BB), is.na(AA), tolerance = 2e-15)
})

#These check out OK.
U4 = test_that("basic pracma::interp2 example gives same outputs large output grid", {

  x <- linspace(-1, 4, 101)
  y <- linspace(-1, 4, 101)
  mgrid <- meshgrid(x, y)
  Z <- mgrid$X^2 + 2*mgrid$Y + 3
  xout <- linspace(-1, 4, 200)
  meshOut <- meshgrid(xout)
  xp <- meshOut$X
  yp <- meshOut$Y
  method <- "linear"
  system.time({AA = interp2(x, y, Z, xp, yp, method)})
  system.time({BB = interp2lin( x, y, Z, xp, yp)})
  expect_equal(sum(AA), sum(BB), tolerance = 2e-15)
})

#These check out OK.
U5 = test_that("basic pracma::interp2 example gives same outputs large input grid", {

  x <- linspace(-1, 4, 1001)
  y <- linspace(-1, 4, 1001)
  mgrid <- meshgrid(x, y)
  Z <- mgrid$X^2 + 2*mgrid$Y + 3
  xout <- linspace(-1, 4, 100)
  meshOut <- meshgrid(xout)
  xp <- meshOut$X
  yp <- meshOut$Y
  method <- "linear"
  system.time({AA = interp2(x, y, Z, xp, yp, method)})
  system.time({BB = interp2lin( x, y, Z, xp, yp)})
  expect_equal(sum(AA), sum(BB), tolerance = 2e-15)
})
