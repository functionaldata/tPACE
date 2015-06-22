# setwd('misc/') 
load('data/InputFormMllwlskInCpp.RData')
if( !exists('Rmullwlsk') ) {
  library(Rcpp)
  sourceCpp('src/Rmullwlsk.cpp')
}
IN = InputFormMllwlskInCpp
 
library(testthat)

# tolerance is relatively large because we cannot control of 2500 * 1e-16 anyway 
# I have already tried using .inverse instead of LLT for the solution and that
# did not make a difference numericallly (small systems anyway)

# These check out OK.
U = test_that("basic Epanetchnikov kernel inputs match MATLAB output for different bandwidths", { 

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = FALSE)
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = FALSE)

  expect_equal(sum(AA), -3.777751915321487e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.540768950285936e+02, tolerance = 1e-15,scale = 1)
  expect_equal(sum(CC),  -3.761635853631063e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
V = test_that("basic rectangular kernel inputs match MATLAB output for different bandwidths", { 

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38), bwCheck = FALSE)
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38), bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38), bwCheck = FALSE)

  expect_equal(sum(AA), -3.288244527254398e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.333882366681741e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.732842331060850e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
H = test_that("basic gaussian kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), bwCheck = FALSE)
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), bwCheck = FALSE)

  expect_equal(sum(AA), -3.689730466900281e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.588111399081811e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.745624473416967e+02, tolerance = 1e-11,scale = 1)
})

# These check out OK.
F = test_that("basic quartic kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38), bwCheck = FALSE)
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38), bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38), bwCheck = FALSE)

  expect_equal(sum(AA), -4.021549928032208e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.472853895316415e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.784289764092692e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
G = test_that("basic gausvar kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38), bwCheck = FALSE)
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38), bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38), bwCheck = FALSE)

  expect_equal(sum(AA), -3.513625436558207e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.480022325255190e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(CC), -3.747722240234172e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
S = test_that("strictly positive window weights inputs match MATLAB output for different bandwidths/kernels", {

  AA = Rmullwlsk( c(4,2), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss', win=seq(1,38)+0, bwCheck = FALSE)
  BB = Rmullwlsk( c(4,2), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect', win=seq(1,38)+0, bwCheck = FALSE)
  CC = Rmullwlsk( c(13.3,23.3), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar', win=seq(1,38)+0, bwCheck = FALSE)
  DD = Rmullwlsk( c(3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan', win=sin(seq(1,38))+3, bwCheck = FALSE)
  EE = Rmullwlsk( c(10.3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar', win=sin(seq(1,38))+3, bwCheck = FALSE)

  expect_equal(sum(AA), -4.337814512732364e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.070237338734315e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(CC), -4.782099399727079e+02, tolerance = 1e-13,scale = 1)
  expect_equal(sum(DD), -4.149192163656717e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(EE), -3.573361421232184e+02, tolerance = 1e-12,scale = 1)
})


# These check out OK.
T = test_that("incoherent kernel_types fall back to Epanechnikov kernels and give correct warnings", {

  AA = Rmullwlsk( c(3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan', win=sin(seq(1,38))+3, bwCheck = FALSE)
  aa = Rmullwlsk( c(3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='boom3', win=sin(seq(1,38))+3, bwCheck = FALSE)

  expect_equal(sum(AA), sum(aa), tolerance = 1e-15,scale = 1)
  expect_warning(  Rmullwlsk(5.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epaffn',win=rep(1,38), bwCheck = FALSE), "Kernel_type argument was not set correctly; Epanechnikov kernel used.")

})



Y = test_that("Small bandwidths give correct error", {

 # expect_error(  Rmullwlsk(0.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = FALSE), "No enough points in local window, please increase bandwidth.")

expect_equal(as.numeric( Rmullwlsk(0.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = TRUE)), 0)
expect_equal(as.numeric( Rmullwlsk(2.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = TRUE)), 1)

expect_equal(as.numeric( Rmullwlsk(0.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = 1)), 0)
expect_equal(as.numeric( Rmullwlsk(2.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = 1)), 1)

expect_equal(as.numeric( Rmullwlsk(0.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = 1)), 0)
expect_equal(as.numeric( Rmullwlsk(0.2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = TRUE)), 0)

})

