devtools::load_all()
library(testthat)
 
test_that('cumtrapzRcpp works on a trivial example', {
  x = c(0,2)
  y = c(0,2)  
  expect_equal( cumtrapzRcpp(x,y), c(0,2) ) 
})

test_that('trapzRcpp works on a nearly trivial example', { 
  x = seq(0,4, length.out=100)
  y = x + sin(x);
  expect_equal(sum( cumtrapzRcpp(x,y)) ,  3.865524746134088e+02) 
})
