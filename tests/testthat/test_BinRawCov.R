devtools::load_all()
library(testthat)

# GMeanAndGCV

set.seed(1)
pts <- c(0, 1, 3:100) / 100
regGrid <- seq(0, 1, by=0.1)
samp3 <- wiener(200, pts, sparsify=2:7)
p0 <- SetOptions(samp3$yList, samp3$tList, optns=list(regular='Sparse', error=TRUE, kernel='epan'))
mu3 <- rep(0, length(pts))
rcov3 <- GetRawCov(samp3$yList, samp3$tList, pts, mu3, p0$regular, error=p0$error)


brcov3 <- BinRawCov(rcov3)
