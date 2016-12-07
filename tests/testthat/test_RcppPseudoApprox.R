# devtools::load_all()
library(testthat)

test_that('RcppPseudoApprox works on a nearly trivial example', {
  set.seed(111)
  z = runif(44);
  expect_equal( RcppPseudoApprox(X = c(0,1), Y = c(0,2) , X_target = z), 2*z, tolerance = 1e-7)
  
}) 

test_that('RcppPseudoApprox errs on obviously wrong data.', { 
  
  expect_equal( RcppPseudoApprox(X = c(0,1), Y = c(0,2,4) , X_target = 0), 2*z, tolerance = 1e-7)
  expect_error( RcppPseudoApprox(X = c(0,1), Y = c(0,2,4) , X_target = 0), 
                 "Problem with unequal vector sizes when doing linear interpolation.")
  
}) 