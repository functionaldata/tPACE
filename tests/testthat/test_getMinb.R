devtools::load_all()
library(testthat)
load('../data/dataForGcvLwlsTest.RData')
# rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output
test_that('2D min bandwidth is similar to Matlab', {
  # expect_equal(getMinb(rcov, sort(unique(unlist(t)))), 4.1427, tolerance=diff(range(unlist(t))) / 1000)
   expect_equal(getMinb(t, sort(unique(unlist(t)))), 4.1427, tolerance=diff(range(unlist(t))) / 1000)
})

test_that('2D min bandwidth for binned and unbinned rcov is the same', {
  # expect_equal(getMinb(BinRawCov(rcov), sort(unique(unlist(t)))), getMinb(rcov, sort(unique(unlist(t)))))
})


set.seed(1)
pts <- seq(0, 1, length=10)
samp2 <- wiener(100, pts, sparsify=2:5)
rcov2 <- GetRawCov(samp2$yList, samp2$tList, pts, rep(0, length(pts)), 'Sparse', FALSE)
getMinb(rcov2, pts)
