cat("\ntests for 'cvlwls1d'")

test_that("basic arguments match MATLAB output ", {

  load('data/dataGeneratedByExampleSeed123.RData')
  a_result = cvlwls1d(x=y,t=t, kernel=kernel, npoly=npoly, nder= nder, dataType='Sparse')
  expect_equal( a_result, 4.172873877723954)

}
)

