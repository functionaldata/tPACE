 cat("\nTests for 'CheckData'")

test_that("basic valid lists arguments do not return any errors ", { 
  CheckData(y = list(c(1,2,3), c(1,2)), t = list(c(1,2,3), c(1,2)))
  CheckData(t = list(c(1,2,3), c(1,2)), y = list(runif(3), runif(2)))
})


#test_that("basic invalid lists arguments do return errors ", { 
#  expect_equal(CheckData( list(c(1,2,3), c(1,2)),list(c(1,2,NA), c(1,2))), TRUE) 
#  expect_equal(CheckData( list(c(1,2,3), c(1,NA)),list(runif(3), runif(2))), TRUE) 
#})

test_that("basic invalid nolists arguments do return errors ", { 
#  expect_equal(CheckData( runif(4), list(c(1,2,NA), c(1,2)) ), TRUE) # We handle NA now
#  expect_equal(CheckData( list(c(1,2,3), c(1,NA)), runif(3) ), TRUE) # We handle NA now
  expect_error(CheckData( matrix(1:6,2,3), list(c(1,2,3), c(1,2)) ), 'y should be list') 
  expect_error(CheckData( list(c(1,2,3), c(1,2)), matrix(1:6,2,3) ), 't should be list') 
})
# cat("Done")
