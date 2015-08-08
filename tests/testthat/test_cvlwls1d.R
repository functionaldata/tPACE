cat("\ntests for 'cvlwls1d'")

test_that("basic arguments match MATLAB output ", {

  try(silent=TRUE, load('data/dataGeneratedByExampleSeed123.RData'))
  try(silent=TRUE, load('tPACE/data/dataGeneratedByExampleSeed123.RData'))

  a_result = cvlwls1d(yy=y, t=t, kernel='epan', npoly=1, nder=0, dataType='Sparse')
  expect_equal( a_result, 4.172873877723954)

}
)

