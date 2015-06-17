# setwd('misc/') 

library(Rcpp)
sourceCpp('src/RrotatedMullwlsk.cpp')
load('data/InputForRotatedMllwlskInCpp.RData')

library(testthat)

# tolerance is relatively large because we cannot control of 2500 * 1e-16 anyway 
# I have already tried using .inverse instead of LLT for the solution and that
# did not make a difference numericallly (small systems anyway)

U = test_that("basic Epanetchnikov kernel inputs match MATLAB output for different bandwidths", { 
  
  AA = Rrotatedmullwlsk(bw =IN$bw, tPairs=IN$tPairs, cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= IN$kernel)
  BB = Rrotatedmullwlsk(bw = c(3,4), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= IN$kernel)
  CC = Rrotatedmullwlsk(bw = c(13,23.3), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= IN$kernel)

  expect_equal(sum(AA), -1.887451898050793, tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -3.264859562745997, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -5.650324984396344, tolerance = 1e-13,scale = 1)
}) 

V = test_that("basic rectangular kernel inputs match MATLAB output for different bandwidths", { 

  AA = Rrotatedmullwlsk(bw =IN$bw, tPairs=IN$tPairs, cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'rect')
  BB = Rrotatedmullwlsk(bw = c(3,4), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'rect')
  CC = Rrotatedmullwlsk(bw = c(13,23.3), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'rect')

  expect_equal(sum(AA),  0.408929466844517, tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -1.803538175275243, tolerance = 1e-13,scale = 1)
  expect_equal(sum(CC), -5.866207150638594, tolerance = 1e-13,scale = 1)

 })

H = test_that("basic gaussian kernel inputs match MATLAB output for different bandwidths", {

  AA = Rrotatedmullwlsk(bw =IN$bw, tPairs=IN$tPairs, cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gauss')
  BB = Rrotatedmullwlsk(bw = c(3,4), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gauss')
  CC = Rrotatedmullwlsk(bw = c(13,23.3), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gauss')

  expect_equal(sum(AA), -4.197686977022681, tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -4.134314374205185, tolerance = 1e-14,scale = 1)
  expect_equal(sum(CC), -5.767647736432502, tolerance = 1e-13,scale = 1)

})

F = test_that("basic quartic kernel inputs match MATLAB output for different bandwidths", {

  AA = Rrotatedmullwlsk(bw =IN$bw, tPairs=IN$tPairs, cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'quar')
  BB = Rrotatedmullwlsk(bw = c(3,4), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'quar')
  CC = Rrotatedmullwlsk(bw = c(13,23.3), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'quar')

  expect_equal(sum(AA), -3.753442160580053,tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -4.970567279909929, tolerance = 1e-13,scale = 1)
  expect_equal(sum(CC), -5.443792883622939, tolerance = 1e-13,scale = 1)

 })




# These check out OK.
G = test_that("basic gausvar kernel inputs match MATLAB output for different bandwidths", {
 
  AA = Rrotatedmullwlsk(bw =IN$bw, tPairs=IN$tPairs, cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gausvar')
  BB = Rrotatedmullwlsk(bw = c(3,4), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gausvar')
  CC = Rrotatedmullwlsk(bw = c(13,23.3), tPairs=(IN$tPairs), cxxn= IN$cxxn, win= IN$win, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gausvar')

  expect_equal(sum(AA), -9.228691155965564, tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -3.594812776733668, tolerance = 1e-13,scale = 1)
  expect_equal(sum(CC), -5.718225024334538, tolerance = 1e-13,scale = 1)

})

# These check out OK.
S = test_that("strictly positive window weights inputs match MATLAB output for different bandwidths/kernels", {

  AA = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= seq(1,38)+0, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gausvar')
  BB = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= seq(1,38)+0, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'gauss')
  CC = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= seq(1,38)+0, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'rect')
  DD = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'epan')
  EE = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'quar')

  expect_equal(sum(AA), -4.924560108566402, tolerance = 1e-13,scale = 1)
  expect_equal(sum(BB), -6.577000474589042, tolerance = 1e-13,scale = 1)
  expect_equal(sum(CC), -1.791956888763226, tolerance = 1e-13,scale = 1)
  expect_equal(sum(DD), -3.614424355861832, tolerance = 1e-13,scale = 1)
  expect_equal(sum(EE), -5.450343839504677, tolerance = 1e-13,scale = 1)
 
})


# These check out OK.
T = test_that("incoherent kernel_types fall back to Epanechnikov kernels and give the proper warning msg.", {

  DD = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'epan')
  dd = Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'boom3')

  expect_equal(sum(DD), sum(dd), tolerance = 1e-15, scale= 1)
  expect_warning(  Rrotatedmullwlsk(bw =c(3,4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'boom3'), "Kernel_type argument was not set correctly; Epanechnikov kernel used.")

})


Y = test_that("Small bandwidths give correct error", {

  expect_error( Rrotatedmullwlsk(bw =c(0.3,0.4), tPairs=IN$tPairs, cxxn= IN$cxxn, win= sin(seq(1,38))+3, xygrid=IN$xygrid, npoly=IN$npoly, kernel_type= 'rect'), "No enough points in local window, please increase bandwidth.")

})

