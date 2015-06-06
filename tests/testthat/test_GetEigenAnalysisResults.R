devtools::load_all()
library(testthat)


set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.01)
samp3 <- wiener(50, pts, sparsify=length(pts))
mu3 <- rep(0, length(pts))

# without error
p0 <- SetOptions(samp3$yList, samp3$tList, list(maxK=50, FVEthreshold=1, dataType='Sparse', error=FALSE, kernel='epan'))
noErrBin <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p0, useBins=TRUE)
tmp <- GetEigenAnalysisResults(noErrBin$smoothCov, regGrid, p0)

debug(GetSmoothedCovarSurface)
debug(GetEigenAnalysisResults)

# with error
p1 <- SetOptions(samp3$yList, samp3$tList, CreateOptions(dataType='Sparse', error=TRUE, kernel='epan'))
Err <- GetSmoothedCovarSurface(samp3$yList, samp3$tList, mu3, pts, regGrid, p1, useBins=FALSE)

# TEst integrate to one.
