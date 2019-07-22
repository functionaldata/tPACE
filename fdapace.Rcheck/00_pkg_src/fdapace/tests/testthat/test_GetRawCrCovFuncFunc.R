# devtools::load_all()
library(testthat)
#source('GetRawCrCovFuncFunc.R')

# These check out OK.
test_that("basic R output matche MATLAB output for different means", {

  # simplest case
  AA = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3)), Lt1 = list(c(1,2,3)), Ymu1 = c(1,4,5), 
                            Ly2 = list(c(1,5)),   Lt2 = list(c(1,2)),   Ymu2 = c(6,5))
  # t_x = {[1 2 3]}; x= {[5 0 3]}; y ={[1 5]}; t_y={[1 2]}; mu_x = [1,4, 5]; mu_y = [6 5]; isYFun=1; regular=0;

  # simplest case different means
  BB = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3)), Lt1 = list(c(1,2,3)), Ymu1 = c(1,4,15), 
                            Ly2 = list(c(1,5)),   Lt2 = list(c(1,2)),   Ymu2 = c(6,15))
  # t_x = {[1 2 3]}; x= {[5 0 3]}; y ={[1 5]}; t_y={[1 2]}; mu_x = [1,4, 15]; mu_y = [6 15]; isYFun=1; regular=0;

  # simple case with single value elements
  CC = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3),c(4)), Lt1 = list(c(1,2,3),c(3)), Ymu1 = c(1,4,15),
                            Ly2 = list(c(1,5),c(5)),   Lt2 = list(c(1,2),c(1)),   Ymu2= c(6,15))
  # t_x = {[1 2 3],[3]}; x= {[5 0 3],[4]}; y ={[1 5],[5]}; t_y={[1 2],[1]}; mu_x = [1,4, 15]; mu_y = [6 15]; isYFun=1; regular=0;

  # simple case with just two vectors each
  DD = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3),c(1,4,8)), Lt1 = list(c(1,2,3),c(1,3,8)), Ymu1 = c(1,4,15,16),
                            Ly2 = list(c(1,5),c(5)),   Lt2 = list(c(1,2),c(1)),   Ymu2= c(6,15))
  # t_x = {[1 2 3],[1 3 8]}; x= {[5 0 3],[1 4 8]}; y ={[1 5],[5]}; t_y={[1 2],[1]}; mu_x = [1,4, 15, 16]; mu_y = [6 15]; isYFun=1; regular=0;
  
  # simple case with just two vectors but extended grid
  EE = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3),c(1,4,8)), Lt1 = list(c(1,2,3),c(1,3,8)), Ymu1 = c(1,4,15,16),
                            Ly2 = list(c(1,5),c(5)),   Lt2 = list(c(1,2),c(16)),  Ymu2= c(6,15,15))
  # t_x = {[1 2 3],[1 3 8]}; x= {[5 0 3],[1 4 8]}; y ={[1 5],[5]}; t_y={[1 2],[16]}; mu_x = [1,4, 15, 16]; mu_y = [6 15 15]; isYFun=1; regular=0; 

  # simple case with three vectors each many duplicate elements
  FF = GetRawCrCovFuncFunc( Ly1 = list(c(5,0,3),c(1,4,8), c(1,2)), Lt1 = list(c(1,2,3),c(1,3,8), c(1,2)), Ymu1 = c(1,4,15,16),
                            Ly2 = list(c(1,5),c(1:5), c(11, 31)),  Lt2 = list(c(1,2),c(16:20), c(1,2)),  Ymu2= c(6,15,15,0,0,0,0))
  # t_x = {[1 2 3],[1 3 8],[1 2]}; x= {[5 0 3],[1 4 8], [1 2]}; mu_x = [1,4, 15, 16]; t_y = {[1 2],  [16:20] [1 2]}; y ={[1 5],  [1:5], [11 31]}; mu_y = [6 15 15 0 0 0 0 ]; isYFun=1; regular=0;

  expect_equal( as.vector(AA$tpairn), c(1, 1, 2, 2, 3, 3, 1, 2, 1, 2, 1, 2) )
  expect_equal( AA$rawCCov, c(-20, 0, 20, 0, 10, 0))

  expect_equal( as.vector(BB$tpairn), c(1, 1, 2, 2, 3, 3, 1, 2, 1, 2, 1, 2) )
  expect_equal( BB$rawCCov, c(-20, -40, 20, 40, 60, 120))

  expect_equal( as.vector(CC$tpairn), c(1, 1, 2, 2, 3, 3, 3, 1, 2, 1, 2, 1, 2, 1) )
  expect_equal( CC$rawCCov, c(-20, -40, 20, 40, 60, 120, 11))

  expect_equal( as.vector(DD$tpairn), c(1, 1, 2, 2, 3, 3, 1, 3, 8, 1, 2, 1, 2, 1, 2, 1, 1, 1) )
  expect_equal( DD$rawCCov, c(-20, -40, 20, 40, 60, 120, 0, 11, 8))

  expect_equal( as.vector(EE$tpairn), c( 1, 1, 2, 2, 3, 3, 1, 3, 8, 1, 2, 1, 2, 1, 2, 16, 16, 16 ) )
  expect_equal( EE$rawCCov, c(-20, -40, 20, 40, 60, 120, 0, 110, 80))

  expect_equal( as.vector(FF$tpairn), c(1, 1, 2, 2, 3, 3, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 8, 8, 8, 8, 8, 1, 1, 2, 2, 1, 2, 1, 2, 1, 2, 16, 17, 18, 19, 20, 16, 17, 18, 19, 20, 16, 17, 18, 19, 20, 1, 2, 1, 2 ) )
  expect_equal( FF$rawCCov, c( -20, -40, 20, 40, 60, 120, 0, 0, 0, 0, 0, 154, -22, -33, -44, -55, 112, -16, -24, -32, -40, 0, 0, -10, -32))
})

test_that('Binned rawCC is the same as unbinned', {
  Ly1 <- list(c(5,0,3),c(1,4,8), c(1,2))
  Lt1 <- list(c(1,2,3),c(1,3,8), c(1,2))
  Ymu1 <- c(1,4,15,16)
  Ly2 <- list(c(1,5),c(1:5), c(11, 31))
  Lt2 <- list(c(1,2),c(16:20), c(1,2))
  Ymu2 <- c(6,15,15,0,0,0,0)
  FF <- GetRawCrCovFuncFunc(Ly1, Lt1, Ymu1, Ly2, Lt2, Ymu2)
  FFbin <- BinRawCov(FF)

  expect_equal(weighted.mean(FFbin$meanVals, FFbin$count), mean(FF$rawCCov))
  expect_equal(sort(unique(FFbin$tPairs[, 1])), sort(unique(FF$tpairn[, 1])))
  expect_equal(sort(unique(FFbin$tPairs[, 2])), sort(unique(FF$tpairn[, 2])))
})
