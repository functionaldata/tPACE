 cat("\nTests for 'CheckData'")

test_that("basic valid lists arguments do not return any errors ", { 
  expect_equal(CheckData( list(c(1,2,3), c(1,2)),list(c(1,2,3), c(1,2))), FALSE) 
  expect_equal(CheckData( list(c(1,2,3), c(1,2)),list(runif(3), runif(2))), FALSE) 
})


test_that("basic invalid lists arguments do return errors ", { 
  expect_equal(CheckData( list(c(1,2,3), c(1,2)),list(c(1,2,NA), c(1,2))), TRUE) 
  expect_equal(CheckData( list(c(1,2,3), c(1,NA)),list(runif(3), runif(2))), TRUE) 
})

test_that("basic invalid nolists arguments do return errors ", { 
  expect_equal(CheckData( runif(4), list(c(1,2,NA), c(1,2)) ), TRUE) 
  expect_equal(CheckData( list(c(1,2,3), c(1,NA)), runif(3) ), TRUE) 
  expect_equal(CheckData( matrix(1:6,2,3), list(c(1,2,3), c(1,2)) ), TRUE) 
  expect_equal(CheckData( list(c(1,2,3), c(1,2)), matrix(1:6,2,3) ), TRUE) 
})
# cat("Done")
