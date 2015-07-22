devtools::load_all()
library(testthat)
 
test_that('trapzRcpp works on a trivial example', {
  x = c(0,2)
  y = c(0,2)  
  expect_equal( trapzRcpp(x,y), 2 ) 
})

test_that('trapzRcpp works on a nearly trivial example', { 
  x = seq(0,4, length.out=100)
  y = x + sin(x);
  expect_equal(trapzRcpp(x,y),  9.653418652171286) 
})
