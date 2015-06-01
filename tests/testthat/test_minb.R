 cat("\nTests for 'minb'")

test_that("basic vector arguments do not return any errors ", { 
  expect_equal( minb(c(1,2,3,4.1),  -1),  NaN) 
  expect_equal( minb(c(1,2,3,4.1),   1),  0.55)
  expect_equal( minb(c(11,2,3,4.1), -1),  NaN)            
  expect_equal( minb(c(11,2,3,4.1),  2),  6.90)
  expect_equal( minb(c(1,2,3,4.1),   6),  NaN)
})

# cat("Done")
