

 # setwd('misc/')
 load('InputFormMllwlskInCpp.RData') 
 library(Rcpp)
 sourceCpp('Rmullwlsk.cpp')
 IN = InputFormMllwlskInCpp
 
library(testthat)

# tolerance is relatively large maybe I still have a bug
# I have already tried using .inverse instead of LLT for the solution and that
# did not make a difference numericallly (small systems anyway)

# These check out OK.
U = test_that(" basic Epanetchnikov inputs match MATLAB output for different bandwidths", { 

AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))
BB = Rmullwlsk( c(5,3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))
CC = Rmullwlsk( c(13.3,23.3),t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38))

  expect_equal(sum(AA ) ,-3.777751915321487e+02 , tolerance = 1e-12,scale = 1)
  expect_equal(sum(BB ) , -3.540768950285936e+02, tolerance = 1e-15,scale = 1)
  expect_equal(sum(CC ) ,  -3.761635853631063e+02, tolerance = 1e-12,scale = 1)
})

 
