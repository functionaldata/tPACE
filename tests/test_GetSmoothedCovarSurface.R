devtools::load_all()
library(testthat)

# GMeanAndGCV

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
samp3 <- wiener(50, pts, sparsify=length(pts))
mu3 <- rep(0, length(pts))

# without error
p0 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(regular='Sparse', error=FALSE, kernel='epan'))
noErr <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p0, useBins=FALSE)

# with error
p1 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(regular='Sparse', error=TRUE, kernel='epan'))
Err <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p1, useBins=FALSE)


# cross-sectional
tmp1 <- do.call(rbind, samp3$yList)
sum((diag(tmp$smoothCov) - seq(0, 1, by=0.1))^2)
sum((diag(cov(tmp1))[seq(1, 21, by=2)] - seq(0, 1, by=0.1))^2)

# GCV
p2 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(bwxcov_gcv='GCV', regular='Sparse', error=FALSE, kernel='epan'))
tmp2 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p2, useBins=FALSE)
sum((diag(tmp2$smoothCov) - seq(0, 1, by=0.1))^2)

# CV
p3 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(bwxcov_gcv='CV', regular='Sparse', error=FALSE, kernel='epan'))
system.time(tmp3 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p3, useBins=FALSE))
sum((diag(tmp3$smoothCov) - seq(0, 1, by=0.1))^2)

# unit tests: test the interface.