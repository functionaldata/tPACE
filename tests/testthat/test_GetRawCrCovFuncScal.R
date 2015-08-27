library(testthat)
source('GetRawCrCovFuncScal.R')

# These check out OK.
test_that("basic R output matche MATLAB output for different means", {

  # simplest case
  AA = GetRawCrCovFuncScal(list(c(5,5,5),c(2,2,2)), list(c(1,2,3), c(1,3,5)), c(1,2,3,4),Z =c(4,5), Zmu = 4)
  # t_x = {[1,2,3], [1,3,5]} ; x= {[5 5 5], [2 2 2]}; y =[ 4 5]; t_y=[]; mu_x = [1,2,3,4]; mu_y = 4; isYFun=0; regular=0;
  
  # simplest case different E[Y]
  BB = GetRawCrCovFuncScal(list(c(5,5,5),c(2,2,2)), list(c(1,2,3), c(1,3,5)), c(1,2,3,4),Z =c(4,5), Zmu = 6)
  # t_x = {[1,2,3], [1,3,5]} ; x= {[5 5 5], [2 2 2]}; y =[ 4 5]; t_y=[]; mu_x = [1,2,3,4]; mu_y = 6; isYFun=0; regular=0;
  
  # simple case more readings per sample
  CC = GetRawCrCovFuncScal(list(c(5,5,5,0),c(2,2,2,0)), list(c(1,2,3,8), c(1,3,5,8)), c(1,2,3,4,1),Z =c(4,5), Zmu =6)
  # t_x = {[1,2,3,8], [1,3,5,8]} ; x= {[5 5 5 0], [2 2 2 0]}; y =[ 4 5]; t_y=[]; mu_x = [1,2,3,4,1]; mu_y = 6; isYFun=0; regular=0;
  
  # simple case more three curves
  DD = GetRawCrCovFuncScal(list(c(5,5,5,0),c(2,2,2,0),c(1,2,5)), list(c(1,2,3,8), c(1,3,5,8), c(1,2,5)), c(1,2,3,4,1),Z =c(4,5,0), Zmu = 0) 
  # t_x = {[1,2,3,8],[1,3,5,8],[1,2,5]} ; x= {[5 5 5 0],[2 2 2 0],[1 2 5]}; y =[ 4 5 0]; t_y=[]; mu_x = [1,2,3,4,1]; mu_y = 0; isYFun=0; regular=0;

  # simple case readings with single measurement
  EE = GetRawCrCovFuncScal(list(c(5,5,5,0),c(2)), list(c(1,2,3,8), c(5)), c(1,2,3,4,1),Z =c(4,5), Zmu = 6)
  # t_x = {[1,2,3,8], [5]} ; x= {[5 5 5 0], [2]}; y =[ 4 5]; t_y=[]; mu_x = [1,2,3,4,1]; mu_y = 6; isYFun=0; regular=0;

  expect_equal( AA$tpairn, c(1,2,3,1,3,5))
  expect_equal( BB$tpairn, c(1,2,3,1,3,5))
  expect_equal( CC$tpairn, c(1,2,3,8,1,3,5,8))
  expect_equal( DD$tpairn, c(1,2,3,8,1,3,5,8,1,2,5)) 
  expect_equal( EE$tpairn, c(1,2,3,8,5))

  expect_equal( AA$raw_ccov, c(0, 0, 0, 1,-1,-2))
  expect_equal( BB$raw_ccov, c(-8,-6,-4,-1,1,2))
  expect_equal( CC$raw_ccov, c(-8,-6,-4,2,-1,1,2,1))
  expect_equal( DD$raw_ccov, c(16,12,8,-4, 5,-5,-10,-5,0,0,0))
  expect_equal( EE$raw_ccov, c(-8,-6,-4,2,2))

})
