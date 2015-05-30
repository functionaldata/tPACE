devtools::load_all()
library(testthat)

p0 <- SetOptions(regular='Sparse', error=TRUE, kernel='epan')

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)
samp3 <- wiener(2000, pts, sparsify=20)
mu3 <- rep(0, length(pts))
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, p0$regular, error=p0$error)
tmp <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p0, useBins=FALSE)
diag(tmp$smoothCov)