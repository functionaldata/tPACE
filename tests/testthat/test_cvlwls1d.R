cat("\ntests for 'CVLwls1D'")

test_that("basic arguments match MATLAB output ", {

  try(silent=TRUE, load('data/dataGeneratedByExampleSeed123.RData'))
  try(silent=TRUE, system.file('testdata', 'dataGeneratedByExampleSeed123.RData', package='fdapace'))

  a_result = CVLwls1D(y, t=t, kernel='epan', npoly=1, nder=0, dataType='Sparse')
  expect_equal( a_result, 4.172873877723954, tol = 0.6) # High tolerance because we have different implementation

}
)

