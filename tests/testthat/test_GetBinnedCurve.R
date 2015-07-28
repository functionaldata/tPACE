devtools::load_all()
library(testthat)

test_that('GetBinnedCurve() works on trivial examples', {
  x = 1:100
  y = 2*x
  A = GetBinnedCurve(x,y,M=10)
  expect_equal(  sum(diff(A$midpoint) ) , 89.100000000000009 )
  expect_equal(  std(A$newy), 60.553007081949829)
  B = GetBinnedCurve(x,y,M=33)
  expect_equal(  sum(diff(B$midpoint) ) , 96 )
  expect_equal(  std(B$newy), 58.069185486196574)
})

test_that('GetBinnedCurve() works on a nearly trivial example', {
  x = seq(0,4, length.out=100)
  y = x + sin(x);
  A = GetBinnedCurve(x,y, 32, TRUE, TRUE, c(1,2.5))
  expect_equal(  sum(diff(A$midpoint) ), 1.453125000000000 )
  expect_equal( std(A$newy), 0.374105692271909  )
})

