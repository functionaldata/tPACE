 cat("\nTests for 'SetOptions'")

 optns = CreateOptions()
 
 test_that("SetOptions test: optns$method = 'IN' for Dense case, 'CE' for Sparse case, 'CE' for other cases with warning", { 
  expect_equal(SetOptions(list(c(1,3,5), c(2,4)),list(c(1,3,5), c(2,4)), optns)$method, 'CE') 
  expect_equal(SetOptions(list(c(1,2,3), c(1,2,3)),list(c(1,2,3), c(1,2,3)), optns)$method, 'IN') 
  expect_equal(SetOptions(list(c(1,2,3,4,5), c(1,2,3,4)),list(c(1,2,3,4,5), c(1,2,3,4)), optns)$method, 'CE') 
})
