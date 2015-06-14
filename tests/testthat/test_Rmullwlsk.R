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

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))

  expect_equal(sum(AA), -3.777751915321487e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.540768950285936e+02, tolerance = 1e-15,scale = 1)
  expect_equal(sum(CC),  -3.761635853631063e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
V = test_that("basic rectangular kernel inputs match MATLAB output for different bandwidths", { 

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38))
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38))
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38))

  expect_equal(sum(AA), -3.288244527254398e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.333882366681741e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.732842331060850e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
H = test_that("basic gaussian kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38))
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38))
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38))

  expect_equal(sum(AA), -3.689730466900281e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.588111399081811e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.745624473416967e+02, tolerance = 1e-11,scale = 1)
})

# These check out OK.
F = test_that("basic quartic kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38))
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38))
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar',win=rep(1,38))

  expect_equal(sum(AA), -4.021549928032208e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.472853895316415e+02, tolerance = 1e-11,scale = 1)
  expect_equal(sum(CC), -3.784289764092692e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
G = test_that("basic gausvar kernel inputs match MATLAB output for different bandwidths", {

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38))
  BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38))
  CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar',win=rep(1,38))

  expect_equal(sum(AA), -3.513625436558207e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.480022325255190e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(CC), -3.747722240234172e+02, tolerance = 1e-12,scale = 1)
})

# These check out OK.
S = test_that("strictly positive window weights inputs match MATLAB output for different bandwidths/kernels", {

  AA = Rmullwlsk( c(4,2), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss', win=seq(1,38)+0)
  BB = Rmullwlsk( c(4,2), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect', win=seq(1,38)+0)
  CC = Rmullwlsk( c(13.3,23.3), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gausvar', win=seq(1,38)+0)
  DD = Rmullwlsk( c(3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan', win=sin(seq(1,38))+3)
  EE = Rmullwlsk( c(10.3,4), t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='quar', win=sin(seq(1,38))+3)

  expect_equal(sum(AA), -4.337814512732364e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB), -3.070237338734315e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(CC), -4.782099399727079e+02, tolerance = 1e-13,scale = 1)
  expect_equal(sum(DD), -4.149192163656717e+02, tolerance = 1e-12,scale = 1)
  expect_equal(sum(EE), -3.573361421232184e+02, tolerance = 1e-12,scale = 1)
})



