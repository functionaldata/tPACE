devtools::load_all()
library(testthat)

p0 <- CreateOptions(regular='Sparse', error=TRUE, kernel='epan')
p1 <- CreateOptions(regular='Dense', error=TRUE, kernel='epan')

set.seed(1)
pts <- seq(0, 1, by=0.05)
regGrid <- seq(0, 1, by=0.1)

sampSparse <- wiener(100, pts, sparsify=20)
p0 <- SetOptions(sampSparse[2], sampSparse[1], p0)
sampDense <- wiener(100, pts, sparsify=21)
p1 <- SetOptions(sampDense[2], sampDense[1], p1)
  
mu3 <- rep(0, length(pts))

## Sparse Case
rcovSparse <- GetRawCov(sampSparse$yList, sampSparse$tList, pts, mu3, p0$regular, error=p0$error)
tmpSparse <- GetSmoothedCovarSurface(sampSparse$yList, sampSparse$tList, mu3, pts, regGrid,
                                     p0, useBins=FALSE)
pccovE_Sparse <- pc_covE(out1 = pts, out21 = seq(0,1,length.out=nrow(tmpSparse$smoothCov)), 
                         bw_xcov = tmpSparse$bwCov, cut = 1, kernel = p0$kernel, rcov = rcovSparse)
pccovE_Sparse$sigma

## Dense Case
rcovDense <- GetRawCov(sampDense$yList, sampDense$tList, pts, mu3, p1$regular, error=p1$error)
tmpDense <- GetSmoothedCovarSurface(sampDense$yList, sampDense$tList, mu3, pts, regGrid,
                                     p1, useBins=FALSE)
pccovE_Dense <- pc_covE(out1 = pts, out21 = seq(0,1,length.out=nrow(tmpDense$smoothCov)), 
                         bw_xcov = tmpDense$bwCov, cut = 1, kernel = p1$kernel, rcov = rcovDense)
pccovE_Dense$sigma
