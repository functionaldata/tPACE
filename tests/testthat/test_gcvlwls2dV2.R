# devtools::load_all()
##options(error=recover)
library(testthat)

try(silent = TRUE, load(system.file('testdata', 'dataForGetRawCov.RData', package='fdapace')))
# try(silent = TRUE, load(system.file('testdata', 'dataForGetRawCov.RData', package='fdapace')))
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) 
r <- range(sort(unlist(t)))
regGrid <- seq(r[1], r[2], length.out=101)
tmp <- GCVLwls2DV2(sort(unlist(t)), regGrid, kern='epan', rcov=rcov, t=t)   ## Why is this guy outside a test?
getGCVscoresV2(tmp$minBW, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid)

test_that('getGCVscoresV2 spits out Inf if bandwidth is too small', {
  expect_equal(suppressWarnings(getGCVscoresV2(2, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid)), Inf)
  expect_message(getGCVscoresV2(2, 'epan', rcov$tPairs, rcov$cxxn, regGrid=regGrid, verbose=TRUE), 'Invalid bandwidth. Try enlarging the window size.\n')
})


# debug(GCVLwls2DV2)
test_that('Binning works for GCV', {
  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  samp3 <- Wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  brcov3 <- BinRawCov(rcov3)
  obsGrid <- sort(unique(unlist(samp3$Lt)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  g1 <- GCVLwls2DV2(obsGrid, regGrid, kern='epan', rcov=rcov3, t=samp3$Lt)
  g2 <- GCVLwls2DV2(obsGrid, regGrid, kern='epan', rcov=brcov3, t=samp3$Lt)
  expect_equal(g1, g2)
})


## Test CV 
test_that('getCVscoresV2 works for binning', {
  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  samp3 <- Wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  brcov3 <- BinRawCov(rcov3)
  obsGrid <- sort(unique(unlist(samp3$Lt)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  partition <- CreateFolds(1:nrow(rcov3$tPairs), k=10)
  bpartition <- CreateFolds(seq_along(brcov3$meanVals), k=10)
  g1 <- getCVscoresV2(partition, 0.35, 'epan', rcov3$tPairs, rcov3$cxxn, regGrid=regGrid)[1]
  g2 <- getCVscoresV2(bpartition, 0.35, 'epan', brcov3$tPairs, brcov3$meanVals, brcov3$count, regGrid=regGrid, RSS=brcov3$RSS)[1]
  expect_equal(g1, g2, tolerance=1e-2)
})

test_that('GCV bandwidth is greater than that of CV', {
  set.seed(2)
  pts <- seq(0, 1, by=0.05)
  samp3 <- Wiener(20, pts, sparsify=2:7)
  rcov3 <- GetRawCov(samp3$Ly, samp3$Lt, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  obsGrid <- sort(unique(unlist(samp3$Lt)))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=101)
  gcvRes <- GCVLwls2DV2(obsGrid, regGrid, kern='epan', rcov=rcov3, t=samp3$Lt)
  cvRes <- GCVLwls2DV2(obsGrid, regGrid, kern='epan', rcov=rcov3, CV=10, t=samp3$Lt)
  expect_gte(gcvRes$h-cvRes$h, -0.01)
})


## HOLE example. Test whether the smoother can handle degenerate cases (in the local window all points lie on a line). 

test_that('GCV will avoid spitting out bandwidth that results in degenerate windows', 
{
  set.seed(2)
  n <- 20
  pts <- seq(0, 1, by=0.1)
  samp4 <- Wiener(20, pts, sparsify=length(pts))
  # for the 10-15th observation retain only one
  for (i in 1:n) {
    retain <- sort(c(1, sample(c(2:3, (length(pts) - 3):length(pts)), 1), 4:(length(pts) - 4)))
    samp4$Ly[[i]] <- samp4$Ly[[i]][retain]
    samp4$Lt[[i]] <- samp4$Lt[[i]][retain]
  }
  
  # CreateDesignPlot(samp4$Lt, pts, TRUE, FALSE, 'samp3')

  rcov4 <- GetRawCov(samp4$Ly, samp4$Lt, pts, rep(0, length(pts)), 'Sparse', error=FALSE)
  g4 <- GCVLwls2DV2(pts, pts, kern='epan', rcov=rcov4, t=samp4$Lt)
  expect_true(g4$minBW > 0.29) # by eyeballing
})

## To matlab
# names(samp3$Lt) <- 1:length(samp3$Lt)
# names(samp3$Ly) <- names(samp3$Lt)
# samp3$Lt <- lapply(samp3$Lt, matrix, nrow=1)
# samp3$Ly <- lapply(samp3$Ly, matrix, nrow=1)
# R.matlab::writeMat('samp3.mat', y=samp3$Ly, t=samp3$Lt)
# GCV values matches matlab, but procedures for optimal GCV BW choice are different.

