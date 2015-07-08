devtools::load_all()
options(error=recover)
library(testthat)

load('../../data/dataForGetRawCov.RData')
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) 
r <- range(sort(unlist(t)))
regGrid <- seq(r[1], r[2], length.out=101)
tmp <- gcvlwls2dV2(sort(unlist(t)), regGrid, kern='epan', rcov=rcov, t=t)
getGCVscoresV2(tmp$minBW, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid)

test_that('getGCVscoresV2 spits out Inf if bandwidth is too small', {
  expect_equal(suppressWarnings(getGCVscoresV2(2, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid)), Inf)
  expect_warning(getGCVscoresV2(2, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid), 'Invalid bandwidth. Try enlarging the window size.')
})


# debug(gcvlwls2dV2)
test_that('Binning works for GCV', {
  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  samp3 <- wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  brcov3 <- BinRawCov(rcov3)
  obsGrid <- sort(unique(unlist(samp3$tList)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  g1 <- gcvlwls2dV2(obsGrid, regGrid, kern='epan', rcov=rcov3, t=samp3$tList)
  g2 <- gcvlwls2dV2(obsGrid, regGrid, kern='epan', rcov=brcov3, t=samp3$tList)
  expect_equal(g1, g2)
})


## Test CV 
test_that('getCVscoresV2 works for binning', {
  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  samp3 <- wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  brcov3 <- BinRawCov(rcov3)
  obsGrid <- sort(unique(unlist(samp3$tList)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  partition <- caret::createFolds(rcov3$cxxn, k=10)
  bpartition <- caret::createFolds(brcov3$meanVals, k=10)
  g1 <- getCVscoresV2(partition, 0.35, 'epan', rcov3$tPairs, rcov3$cxxn, regGrid=regGrid)
  g2 <- getCVscoresV2(bpartition, 0.35, 'epan', brcov3$tPairs, brcov3$meanVals, brcov3$count, regGrid=regGrid, RSS=brcov3$RSS)
  expect_equal(g1, g2, tolerance=1e-2)
})

test_that('GCV is closed to CV', {
  set.seed(2)
  pts <- seq(0, 1, by=0.05)
  samp3 <- wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  obsGrid <- sort(unique(unlist(samp3$tList)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  gcvRes <- gcvlwls2dV2(obsGrid, regGrid, kern='epan', rcov=rcov3, t=samp3$tList)
  cvRes <- gcvlwls2dV2(obsGrid, regGrid, kern='epan', rcov=rcov3, CV='10fold', t=samp3$tList)
  expect_equal(gcvRes$h, cvRes$h, 0.1)
})


## HOLE example. Test whether the smoother can handle degenerate cases (in the local window all points lie on a line). 

test_that('GCV will avoid spitting out bandwidth that results in degenerate windows', 
{
  set.seed(2)
  n <- 20
  pts <- seq(0, 1, by=0.1)
  samp4 <- wiener(20, pts, sparsify=length(pts))
  # for the 10-15th observation retain only one
  for (i in 1:n) {
    retain <- sort(c(1, sample(c(2:3, (length(pts) - 3):length(pts)), 1), 4:(length(pts) - 4)))
    samp4$yList[[i]] <- samp4$yList[[i]][retain]
    samp4$tList[[i]] <- samp4$tList[[i]][retain]
  }
  createDesignPlot(samp4$tList, pts, TRUE, FALSE, 'samp3')

  rcov4 <- GetRawCov(samp4$yList, samp4$tList, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  g4 <- gcvlwls2dV2(pts, pts, kern='epan', rcov=rcov4, t=samp4$tList)
  expect_true(g4$minBW > 0.3) # by eyeballing
})

## To matlab
# names(samp3$tList) <- 1:length(samp3$tList)
# names(samp3$yList) <- names(samp3$tList)
# samp3$tList <- lapply(samp3$tList, matrix, nrow=1)
# samp3$yList <- lapply(samp3$yList, matrix, nrow=1)
# R.matlab::writeMat('samp3.mat', y=samp3$yList, t=samp3$tList)
# GCV values matches matlab, but procedures for optimal GCV BW choice are different.
