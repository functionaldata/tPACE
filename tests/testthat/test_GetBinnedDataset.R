devtools::load_all()
library(testthat)

test_that('GetBinnedDataset() works on trivial examples', {

  y = list( c(1:1000), c(3:1012))
  t = list( seq(0,1,length.out=1000), seq(0,1,length.out=1010))
  A  = GetBinnedDataset(y,t, optns=SetOptions(y,t,NULL))
 
  expect_equal(   length(A$newt[[2]]), 400)
  expect_equal(   A$newt[[2]][113] ,  0.281250000000000)
  expect_equal( mean( A$newt[[1]]),  mean( A$newt[[2]]) )
  expect_equal(   A$newy[[1]][313] , 7.815000000000000e+02 )
  expect_equal( mean( A$newy[[1]]), 5.005000000000000e+02 )

})
