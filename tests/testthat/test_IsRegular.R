cat("\ntests for 'IsRegular'")

test_that("basic valid lists arguments do not return any errors", { 
  expect_equal(IsRegular(list(c(1:10), c(1:10), c(1:10))), 'Dense') 
  expect_equal(IsRegular(list(c(1,2,3  ), c(1,2,3,4), c(1,2,3,4))), 'DenseWithMV') 
  expect_equal(IsRegular(list(c(1,2   ), c(1,2,3  ), c(1,2,3,4))), 'Sparse') 
}
)

cat("\nexcept for dense but irregular case")
IsRegular(list(c(1,2,3,5),c(1,2,3,5),c(1,2,3,5)))
