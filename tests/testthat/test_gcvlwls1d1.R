 cat("\nTests for 'GCVLwls1D1.R'")

load('data/dataForGcvLwlsTest.RData')

test_that("basic  optimal bandwidth choice for the mean function use GCV method matches MATLAB for Sparse data", {

  A <- GCVLwls1D1(y,t,'epan',1,0,'Sparse') 
  expect_equal( A$bOpt, 2.071354057811459 ) 

  B <- GCVLwls1D1(y,t,'rect',1,0,'Sparse') 
  expect_equal( B$bOpt, 2.238990337557121, 0.04 ) 

 } )

