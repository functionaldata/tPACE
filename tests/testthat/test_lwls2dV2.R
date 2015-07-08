devtools::load_all()
library(testthat)

# as in test_Rmullwlsk.R
load('../../data/InputFormMllwlskInCpp.RData')
IN = InputFormMllwlskInCpp



test_that('lwls2dV2 interface is correct using xout1 and xout2', {
  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE)
  BB = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE)
  CC = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE)
  DD = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE)
  expect_equal(lwls2dV2(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), AA)
  expect_equal(lwls2dV2(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn), BB)
  expect_equal(lwls2dV2(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout1=IN$regGrid, xout2=IN$regGrid), CC)
  expect_equal(lwls2dV2(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn), DD)
})


test_that('lwls2dV2 interface is correct using xout', {
  AA = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), FALSE))
  BB = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='gauss',win=rep(1,38), FALSE))
  CC = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), FALSE))
  DD = diag(Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=sort(unique(IN$tPairs[, 1])), ygrid=sort(unique(IN$tPairs[, 2])), kernel_type='epan',win=rep(1,38), FALSE))
  expect_equal(lwls2dV2(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), AA, tolerance=0.05)
  expect_equal(lwls2dV2(2* IN$bw, kern='gauss', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), BB)
  expect_equal(lwls2dV2(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(IN$regGrid, IN$regGrid)), CC, tolerance=0.05)
  expect_equal(lwls2dV2(2* IN$bw, kern='epan', IN$tPairs, IN$cxxn, xout=cbind(sort(unique(IN$tPairs[, 1])), sort(unique(IN$tPairs[, 2])))), DD)
})

