#options(error=recover)
library(testthat)

test_that('SimpleFolds works', {
  samp <- 1:10
  res1 <- SimpleFolds(samp, 10)
  expect_equal(length(res1), 10)
  expect_equal(sum(sapply(res1, length)), length(samp))
  expect_equal(diff(range(sapply(res1, length))), 0)

  samp <- 20:60
  res2 <- SimpleFolds(samp, 10)
  expect_equal(length(res2), 10)
  expect_equal(sum(sapply(res2, length)), length(samp))
  expect_equal(diff(range(sapply(res2, length))), 1)
})

test_that('CreateFolds works for numeric responses', {
  set.seed(1)
  samp1 <- CreateFolds(rnorm(55), 10) 
  sampVec1 <- sort(do.call(c, samp1))
  expect_true(diff(range(sapply(samp1, length))) <= 5)
  expect_equal(sampVec1, seq_along(sampVec1))
  
  samp2 <- CreateFolds(rnorm(5), 10) 
  sampVec2 <- sort(do.call(c, samp2))
  expect_true(diff(range(sapply(samp2, length))) <= 1)
  expect_equal(sampVec2, seq_along(sampVec2))
})

test_that('CreateFolds works for factor/class responses', {
  set.seed(1)
  nclass <- 2
  tmp <- sample(c(rep(0, 10), rep(1, 11)))
  samp3 <- CreateFolds(tmp, 10) 
  expect_true(diff(range(sapply(samp3, length))) <= nclass)
  expect_true(all(sapply(samp3, function(ind) any(tmp[ind] == 0) && any(tmp[ind] == 1))))

  nclass <- 3
  tmp <- sample(c(rep(1, 10), rep(2, 11), rep(3, 5)))
  tmp <- factor(tmp)
  samp4 <- CreateFolds(tmp, 10) 
  expect_true(diff(range(sapply(samp4, length))) <= nclass)
  expect_true(all(sapply(samp4, function(ind) sum(table(tmp[ind]) != 0) >= 2)))
})

# debug(CreateFolds)
# undebug(CreateFolds)
# debug(SimpleFolds)
# undebug(SimpleFolds)
