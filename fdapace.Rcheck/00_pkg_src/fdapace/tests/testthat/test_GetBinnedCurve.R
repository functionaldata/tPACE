# devtools::load_all()
library(testthat)

test_that('GetBinnedCurve() works on trivial examples', {
  x = 1:100
  y = 2*x
  A = GetBinnedCurve(x,y,M=50)
  expect_equal(  sum(diff(A$midpoint) ) , diff(range(x)) )
  expect_equal(  sd(A$newy), sd(y), tolerance=0.05)
  B = GetBinnedCurve(x,y,M=33)
  expect_equal(  sum(diff(B$midpoint) ) , diff(range(x)))
  expect_equal(  sd(B$newy), sd(y), tolerance=0.05)
})

test_that('GetBinnedCurve() works on a nearly trivial example', {
  x = seq(0,4, length.out=100)
  y = x + sin(x);
  A = GetBinnedCurve(x,y, 32, TRUE, TRUE, c(1,2.5))
  expect_equal(  sum(diff(A$midpoint) ), 1.5, tolerance=0.05 )
  expect_equal( sd(A$newy), 0.38, tolerance=0.05)
})

test_that('GetBinnedCurve() works for large case',{
  x <- 1:2000
  y <- 2 * x
  A = GetBinnedCurve(x, y, 400, TRUE, TRUE, c(1, 1999))
  expect_equal(  A$binWidth, 4.995 )
  expect_lte(  sd(A$count), 0.2)
  expect_equal(  sd(A$newy), sd(y), tolerance=0.03)
  expect_equal(  sd(A$midpoint), sd(x), tolerance=0.03)
})






