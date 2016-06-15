# devtools::load_all()
library(testthat)
load('data/dataForGcvLwlsTest.RData')
# rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output
test_that('2D min bandwidth is similar to Matlab', {
  # expect_equal(GetMinb(rcov, sort(unique(unlist(t)))), 4.1427, tolerance=diff(range(unlist(t))) / 1000)
  
  # We break strict compatibility with MATLAB when we used  quantile( diff(b[ids]), 0.95)/2 instead of max(diff(b[ids])/2)
  # expect_equal(GetMinb(t, sort(unique(unlist(t)))), 4.1427, tolerance=diff(range(unlist(t))) / 1000)
  
  expect_equal( GetMinb(legacyCode = TRUE,t, sort(unique(unlist(t)))), 4.1427, tolerance= 4.1427 * 0.0175)
  expect_equal( as.numeric(GetMinb(t, sort(unique(unlist(t))))), 4.1427, tolerance= 4.1427 * 0.0195 )
})

test_that('2D min bandwidth for binned and unbinned rcov is the same', {
  # expect_equal(GetMinb(BinRawCov(rcov), sort(unique(unlist(t)))), GetMinb(rcov, sort(unique(unlist(t)))))
})


set.seed(1)
pts <- seq(0, 1, length=10)
samp2 <- Wiener(100, pts, sparsify=2:5)
rcov2 <- GetRawCov(samp2$Ly, samp2$Lt, pts, rep(0, length(pts)), 'Sparse', FALSE)
GetMinb(samp2$Lt, pts)
