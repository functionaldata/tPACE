 cat("\nTests for 'gcvlwls1d1.R'")

load('../../data/dataForGcvLwlsTest.RData')

test_that("basic  optimal bandwidth choice for the mean function use GCV method matches MATLAB for Sparse data", {

  A <- gcvlwls1d1(y,t,'epan',1,0,'Sparse') 
  expect_equal( A$bOpt, 2.071354057811459 ) 

  B <- gcvlwls1d1(y,t,'rect',1,0,'Sparse') 
  expect_equal( B$bOpt, 2.238990337557121 ) 

 } )

