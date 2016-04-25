 cat("\nTests for 'Minb'")

test_that("basic vector arguments do not return any errors ", { 
  expect_equal( Minb(c(1,2,3,4.1),  -1),  NaN) 
  expect_equal( Minb(c(1,2,3,4.1),   1),  2* 0.55)
  expect_equal( Minb(c(11,2,3,4.1), -1),  NaN)            
  expect_equal( Minb(c(11,2,3,4.1),  2),  8)
  expect_equal( Minb(c(1,2,3,4.1),   6),  NaN)
})

# cat("Done")
