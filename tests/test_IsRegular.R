 cat("\nTests for 'IsRegular'")

test_that("basic valid lists arguments do not return any errors ", { 
  expect_equal(IsRegular(list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4))), 2) 
  expect_equal(IsRegular(list(c(1,2,3  ), c(1,2,3,4), c(1,2,3,4))), 1) 
  expect_equal(IsRegular(list(c(1,2   ), c(1,2,3  ), c(1,2,3,4))), 0) 
})

# cat("Done")
