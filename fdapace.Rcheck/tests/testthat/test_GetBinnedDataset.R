# devtools::load_all()
library(testthat)

test_that('GetBinnedDataset() works on trivial examples', {

  y = list( c(1:1000), c(3:1012))
  t = list( seq(0,1,length.out=1000), seq(0,1,length.out=1010))
  A  = GetBinnedDataset(y,t, optns=SetOptions(y,t,NULL))
 
  expect_equal(   length(A$newt[[2]]), 400)
  expect_equal(   A$newt[[2]][113] ,  0.2807, tolerance=0.01, scale=1)
  expect_equal( mean( A$newt[[1]]),  mean( A$newt[[2]]) )
  expect_equal(   A$newy[[1]][313] , 782, tolerance=0.01, scale=1)
  expect_equal( mean( A$newy[[1]]), 5.005000000000000e+02 )

})

test_that('GetBinnedDataset() works on binned examples, with the first individual being singleton', {

  y = list(1,  seq(0, 1000, length.out=1000), seq(3, 1012, length.out=1010))
  t = list(0.5, seq(0,1,length.out=1000), seq(0,1,length.out=1010))
  A  = GetBinnedDataset(y,t, optns=SetOptions(y,t, list(numBins=20)))
 
  expect_equal(   length(A$newt[[3]]), 20)
  expect_true(abs(A$newt[[2]][10] - A$newy[[2]][10] / 1000) < 1e-3)
  expect_equal( mean( A$newt[[3]]),  mean( A$newt[[2]]) )

  })
