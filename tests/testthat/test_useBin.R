devtools::load_all()
library(testthat)

test_that('The binned version is exactly the same as the unbinned version.', {
  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  samp3 <- Wiener(300, pts, sparsify=5)
  y <- samp3$yList
  t <- samp3$tList

  resNoBin <- FPCA(y, t, list(dataType='Sparse', useBinnedData='OFF', useBinnedCov=FALSE))
  resBin <- FPCA(y, t, list(dataType='Sparse', useBinnedData='OFF', useBinnedCov=TRUE))

  expect_equal(resNoBin[names(resNoBin) != 'optns'], resBin[names(resBin) != 'optns'])
})
