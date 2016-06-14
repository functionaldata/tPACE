library(testthat)
load('data/InputFormMllwlskInCpp.RData')
IN = InputFormMllwlskInCpp
 
ord <- order(IN$tPairs[, 1])
xin <- IN$tPairs[ord, ]
yin <- IN$cxxn[ord]
win <- rep(1,38)  


# These check out OK.
U = test_that("basic inputs for different kernels match previous inputs.", { 

  AA = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan',win=rep(1,38), bwCheck = FALSE)
  AAu = RmullwlskUniversal( bw = 2* IN$bw, tPairs =t(xin), cxxn=yin, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='epan', win=win, FALSE, TRUE) 
   
  BB = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect',win=rep(1,38), bwCheck = FALSE)
  BBu = RmullwlskUniversal( bw = 2* IN$bw, tPairs =t(xin), cxxn=yin, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='rect', win=win, FALSE, TRUE) 
  
  CC = Rmullwlsk(2* IN$bw,t(IN$tPairs),cxxn=IN$cxxn, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss',win=rep(1,38), bwCheck = FALSE)
  CCu = RmullwlskUniversal( bw = 2* IN$bw, tPairs =t(xin), cxxn=yin, xgrid=IN$regGrid, ygrid=IN$regGrid, kernel_type='gauss', win=win, FALSE, TRUE) 
  
  expect_equal( CC, CCu)
  expect_equal( BB, BBu)
  expect_equal( AA, AAu)



})


