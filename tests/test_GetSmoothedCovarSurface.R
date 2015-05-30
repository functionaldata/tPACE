devtools::load_all()
library(testthat)

# GMeanAndGCV
p0 <- SetOptions(regular='Sparse', error=FALSE, kernel='epan')

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
samp3 <- wiener(200, pts, sparsify=2:7)
mu3 <- rep(0, length(pts))
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, p0$regular, error=p0$error)
tmp <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p0, useBins=FALSE)
# cross-sectional
tmp1 <- do.call(rbind, samp3$yList)
sum((diag(tmp$smoothCov) - seq(0, 1, by=0.1))^2)
sum((diag(cov(tmp1))[seq(1, 21, by=2)] - seq(0, 1, by=0.1))^2)

# GCV
p2 <- SetOptions(bwxcov_gcv='GCV', regular='Sparse', error=FALSE, kernel='epan')
system.time(tmp2 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p2, useBins=FALSE))
sum((diag(tmp2$smoothCov) - seq(0, 1, by=0.1))^2)

# CV
p3 <- SetOptions(bwxcov_gcv='CV', regular='Sparse', error=FALSE, kernel='epan')
system.time(tmp3 <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p3, useBins=FALSE))
sum((diag(tmp3$smoothCov) - seq(0, 1, by=0.1))^2)