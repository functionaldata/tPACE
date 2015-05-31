library(testthat)
load('data/dataForGcvLwlsTest.RData')
rcov <- GetRawCov(y,t, sort(unlist(t)), mu,'Sparse',FALSE) #Matches ML output
getMinb(rcov, sort(unique(unlist(t)))) # close to 4.1427 given by getMinb.m



set.seed(1)
pts <- seq(0, 1, length=10)
samp2 <- wiener(100, pts, sparsify=2:5)
rcov2 <- GetRawCov(samp2$yList, samp2$tList, pts, rep(0, length(pts)), 'Sparse', FALSE)
getMinb(rcov2, pts)
